"""Logique de la section « ProteoCast des ABP » (garde streamlit.py léger).

Fournit :
  * la liste des ABP + UniProt + statut (résultat ProteoCast déposé ou non)
  * les chaînes d'interface d'un ABP (pour la piste verte du paysage mutationnel)
  * les images du dossier de résultats ProteoCast si le dossier complet est déposé

Le rendu du paysage (_render_abp_proteocast, _render_abp_actin_conservation) reste
dans streamlit.py ; ce module ne fait que préparer les données.
"""
import os
import re
import urllib.request
import pandas as pd
import streamlit as st

_MANIFEST = "data/proteocast/abp_inputs/manifest.csv"
_ABP_DIR = "data/proteocast/abp"
_ALL = "data/filtered/filtered_all_data.csv"
_IFACE = "data/filtered/details/3.interface_residues.csv"


def _slug(s):
    return re.sub(r"[^a-zA-Z0-9]+", "_", str(s)).strip("_")[:60]


def _result_dir_or_csv(slug):
    """Renvoie le dossier de résultats complet si présent, sinon None."""
    d = os.path.join(_ABP_DIR, slug)
    return d if os.path.isdir(d) else None


@st.cache_data(show_spinner=False)
def load_status(_mtime):
    """Table ABP : titre, uniprot, slug, statut (fait / à faire), fusion."""
    if not os.path.exists(_MANIFEST):
        return None
    m = pd.read_csv(_MANIFEST).sort_values("abp_title").reset_index(drop=True)

    def done(row):
        p1 = os.path.join(_ABP_DIR, f"{row['slug']}.csv")
        p2 = os.path.join(_ABP_DIR, row["slug"], "4.query_ProteoCast.csv")
        return os.path.exists(p1) or os.path.exists(p2)
    m["fait"] = m.apply(done, axis=1)
    return m


def run_proteocast_job(uniprot, slug, log=None):
    """Soumet un job ProteoCast (Retrieve from AlphaFoldDB), attend, télécharge et
    range le résultat dans data/proteocast/abp/<slug>/. Réutilise le script
    proteocast_submit_abp. `log` = callback(str) optionnel (feedback Streamlit).
    Renvoie (ok: bool, message: str)."""
    import requests
    import proteocast_submit_abp as pcs

    def _say(t):
        if log:
            log(t)
    sess = requests.Session()
    sess.headers.update({"User-Agent": pcs.UA})
    try:
        jid = pcs.submit(sess, uniprot)
    except Exception as e:
        return False, f"Submission failed: {e}"
    if not jid:
        return False, "No job_id returned by the server."
    _say(f"Job {jid} submitted — computing on proteocast.ijm.fr…")
    status, msg = pcs.wait(sess, jid)
    if status != "finished":
        return False, f"Status “{status}”: {msg or 'no message'}"
    _say("Computation finished — downloading the result folder…")
    try:
        ok = pcs.download_extract(sess, jid, slug)
    except Exception as e:
        return False, f"Download failed: {e}"
    return (ok, "ProteoCast fetched." if ok
            else "ZIP downloaded but key file missing.")


@st.cache_data(show_spinner=False)
def abp_interface_chains(title, _mtime):
    """Chaînes d'interface (casse d'interface_residues) contactées par cet ABP,
    quel que soit le sens actin/ABP. Sert à la piste verte du paysage."""
    if not (os.path.exists(_ALL) and os.path.exists(_IFACE)):
        return set()
    df = pd.read_csv(_ALL, low_memory=False)
    df["s1_actine"] = df["s1_actine"].fillna(False).astype(bool)
    df["s2_actine"] = df["s2_actine"].fillna(False).astype(bool)
    chains = set()
    d1 = df[df["s1_actine"] & ~df["s2_actine"] & (df["subunit_2_title"] == title)]
    chains |= set(d1["subunit_2"].dropna().str.lower())
    d2 = df[~df["s1_actine"] & df["s2_actine"] & (df["subunit_1_title"] == title)]
    chains |= set(d2["subunit_1"].dropna().str.lower())
    if not chains:
        return set()
    ir = pd.read_csv(_IFACE, usecols=["chain"])
    # remonter à la casse réelle d'interface_residues
    return {c for c in ir["chain"].dropna().unique() if str(c).lower() in chains}


def result_images(slug):
    """Images du dossier de résultats ProteoCast (paysage, GMM…) si dossier déposé."""
    d = _result_dir_or_csv(slug)
    if d is None:
        return []
    return sorted(os.path.join(d, f) for f in os.listdir(d)
                  if f.lower().endswith((".png", ".jpg", ".jpeg")))


@st.cache_data(show_spinner="Domaines InterPro…")
def fetch_domains(uniprot):
    """Domaines d'un UniProt via l'API InterPro : liste {name, db, acc, spans}.
    spans = liste de (start, end). Garde Pfam / SMART / InterPro (domain/repeat/
    family). None si échec."""
    import json
    if not uniprot:
        return None
    url = (f"https://www.ebi.ac.uk/interpro/api/entry/all/protein/uniprot/"
           f"{uniprot}/?page_size=200")
    try:
        req = urllib.request.Request(url, headers={"Accept": "application/json"})
        with urllib.request.urlopen(req, timeout=30) as r:
            data = json.load(r)
    except Exception:
        return None
    _KEEP_DB = {"pfam", "smart", "interpro"}
    _KEEP_TYPE = {"domain", "repeat", "family"}
    out = []
    for e in data.get("results", []):
        meta = e.get("metadata", {})
        db = meta.get("source_database", "")
        typ = meta.get("type", "")
        if db not in _KEEP_DB or typ not in _KEEP_TYPE:
            continue
        spans = []
        for p in e.get("proteins", []):
            for loc in p.get("entry_protein_locations", []):
                for frag in loc.get("fragments", []):
                    if frag.get("start") and frag.get("end"):
                        spans.append((int(frag["start"]), int(frag["end"])))
        if spans:
            out.append({"name": meta.get("name") or meta.get("accession"),
                        "db": db, "acc": meta.get("accession"),
                        "type": typ, "spans": spans})
    # dédupe par (name, spans) ; priorité Pfam > SMART > InterPro
    _order = {"pfam": 0, "smart": 1, "interpro": 2}
    out.sort(key=lambda d: (_order.get(d["db"], 9), d["name"]))
    seen, uniq = set(), []
    for d in out:
        key = (d["name"], tuple(sorted(d["spans"])))
        if key in seen:
            continue
        seen.add(key)
        uniq.append(d)
    return uniq


# ASA max par résidu (Tien et al. 2013, valeurs théoriques) pour normaliser en RSA
_MAXASA = {"A": 129, "R": 274, "N": 195, "D": 193, "C": 167, "E": 223,
           "Q": 225, "G": 104, "H": 224, "I": 197, "L": 201, "K": 236,
           "M": 224, "F": 240, "P": 159, "S": 155, "T": 172, "W": 285,
           "Y": 263, "V": 174}
_ABP_CHAINS = "data/exports/abp_site_domain/abp_chains"
_ABP_REPS = "data/exports/abp_site_domain/abp_representatives.csv"


@st.cache_data(show_spinner="Computing RSA from the ABP structure…")
def compute_abp_rsa(abp_title, _mtime):
    """RSA par position, calculé sur la VRAIE chaîne ABP de nos structures
    (SASA Shrake-Rupley, biopython), pas le monomère AlphaFold.
    Renvoie {position (numérotation séquence) : RSA fraction} ou {}."""
    import glob
    import re
    if not (os.path.exists(_ABP_REPS) and os.path.isdir(_ABP_CHAINS)):
        return {}
    reps = pd.read_csv(_ABP_REPS)
    row = reps[reps["abp_title"] == abp_title]
    if row.empty:
        return {}
    r = row.iloc[0]
    slug = re.sub(r"[^A-Za-z0-9]+", "_", str(abp_title)).strip("_")[:50]
    cand = glob.glob(os.path.join(_ABP_CHAINS, f"{slug}__*.pdb"))
    if not cand:
        return {}
    try:
        from Bio.PDB import PDBParser
        from Bio.PDB.SASA import ShrakeRupley
        from Bio.Data.IUPACData import protein_letters_3to1
    except ImportError:
        return {}
    st_ = PDBParser(QUIET=True).get_structure("x", cand[0])
    model = next(iter(st_))
    ShrakeRupley().compute(model, level="R")
    out = {}
    for ch in model:
        for res in ch:
            if res.id[0] != " ":
                continue
            aa1 = protein_letters_3to1.get(res.resname.capitalize(), None)
            if aa1 is None or aa1 not in _MAXASA:
                continue
            rsa = res.sasa / _MAXASA[aa1]
            out[int(res.id[1])] = min(rsa, 1.5)
    return out


def _query_seq(slug):
    """Séquence exacte soumise à ProteoCast (1.query.fasta) ou None."""
    qf = os.path.join(_ABP_DIR, slug, "1.query.fasta")
    if not os.path.exists(qf):
        return None
    q = "".join(l.strip() for l in open(qf) if not l.startswith(">"))
    return q or None


def _align_map(src, q):
    """Map {position 1-based sur `src` -> position 1-based sur `q`} par alignement
    global. Indispensable quand la construction PDB de NOS structures ne suit pas
    la numérotation UniProt soumise à ProteoCast (extension N-term, isoforme,
    résidus manquants…). Identité si séquences identiques ; {} si échec."""
    if not src or not q:
        return {}
    if src == q:
        return {i + 1: i + 1 for i in range(len(q))}
    try:
        from Bio.Align import PairwiseAligner, substitution_matrices
    except ImportError:
        return {}
    al = PairwiseAligner()
    al.mode = "global"
    al.substitution_matrix = substitution_matrices.load("BLOSUM62")
    al.open_gap_score = -10
    al.extend_gap_score = -0.5
    try:
        aln = al.align(src, q)[0]
    except Exception:
        return {}
    mp = {}
    for (s0, s1), (q0, q1) in zip(*aln.aligned):
        for k in range(s1 - s0):
            mp[s0 + k + 1] = q0 + k + 1
    return mp


@st.cache_data(show_spinner=False)
def abp_interface_asa_on_query(abp_title, slug, _mtime):
    """Empreinte d'interface de l'ABP DÉJÀ replacée sur la numérotation query
    ProteoCast : {position query : %ASA enfouie max}.

    residue_number_sequence est 1-based sur la s2_sequence de CHAQUE structure ;
    on aligne donc la séquence de chaque chaîne sur la query, puis on remappe.
    C'est ce qui recale la piste « ABP contact actin » sur la heatmap et les
    domaines (ex. Adducin : construction PDB != UniProt query)."""
    q = _query_seq(slug)
    if not (q and os.path.exists(_ALL) and os.path.exists(_IFACE)):
        return {}
    df = pd.read_csv(_ALL, low_memory=False)
    # chaîne (casse d'interface_residues) -> séquence de cette structure
    chain_seq = {}
    for col_c, col_s, col_t in (("subunit_2", "s2_sequence", "subunit_2_title"),
                                ("subunit_1", "s1_sequence", "subunit_1_title")):
        if not {col_c, col_s, col_t} <= set(df.columns):
            continue
        d = df[df[col_t] == abp_title].dropna(subset=[col_c, col_s])
        for c, s in zip(d[col_c], d[col_s]):
            chain_seq.setdefault(str(c).lower(), s)
    if not chain_seq:
        return {}
    ir = pd.read_csv(_IFACE)
    ir["_low"] = ir["chain"].astype(str).str.lower()
    ir = ir[ir["_low"].isin(chain_seq)]
    if ir.empty:
        return {}
    ir["_asa"] = pd.to_numeric(
        ir["buried_ASA_percent"].astype(str).str.replace("%", "", regex=False),
        errors="coerce")
    ir["_pos"] = pd.to_numeric(ir["residue_number_sequence"], errors="coerce")
    ir = ir.dropna(subset=["_asa", "_pos"])
    maps = {c: _align_map(s, q) for c, s in chain_seq.items()}  # 1 align / chaîne
    out = {}
    for c, pos, asa in zip(ir["_low"], ir["_pos"].astype(int), ir["_asa"]):
        qp = maps.get(c, {}).get(pos)
        if qp:
            out[qp] = max(out.get(qp, 0.0), float(asa))
    return out


@st.cache_data(show_spinner="Computing RSA from the ABP structure…")
def abp_rsa_on_query(abp_title, slug, _mtime):
    """RSA de l'ABP replacé sur la numérotation query ProteoCast :
    {position query : RSA}. Le RSA est calculé sur la VRAIE chaîne ABP de nos
    structures (SASA Shrake-Rupley), puis la séquence extraite du PDB est alignée
    sur la query pour recaler les positions (comme l'empreinte)."""
    import glob
    q = _query_seq(slug)
    if not (q and os.path.exists(_ABP_REPS) and os.path.isdir(_ABP_CHAINS)):
        return {}
    reps = pd.read_csv(_ABP_REPS)
    if reps[reps["abp_title"] == abp_title].empty:
        return {}
    slug50 = re.sub(r"[^A-Za-z0-9]+", "_", str(abp_title)).strip("_")[:50]
    cand = glob.glob(os.path.join(_ABP_CHAINS, f"{slug50}__*.pdb"))
    if not cand:
        return {}
    try:
        from Bio.PDB import PDBParser
        from Bio.PDB.SASA import ShrakeRupley
        from Bio.Data.IUPACData import protein_letters_3to1
    except ImportError:
        return {}
    st_ = PDBParser(QUIET=True).get_structure("x", cand[0])
    model = next(iter(st_))
    ShrakeRupley().compute(model, level="R")
    seq, rsa_by_rank = [], {}      # séquence PDB dans l'ordre + RSA par rang
    for ch in model:
        for res in ch:
            if res.id[0] != " ":
                continue
            aa1 = protein_letters_3to1.get(res.resname.capitalize())
            if aa1 is None or aa1 not in _MAXASA:
                continue
            rsa_by_rank[len(seq) + 1] = min(res.sasa / _MAXASA[aa1], 1.5)
            seq.append(aa1)
        break                       # 1re chaîne protéique
    mp = _align_map("".join(seq), q)
    return {mp[r]: v for r, v in rsa_by_rank.items() if r in mp}


def _proteocast_csv(slug):
    for c in (os.path.join(_ABP_DIR, f"{slug}.csv"),
              os.path.join(_ABP_DIR, slug, "4.query_ProteoCast.csv")):
        if os.path.exists(c):
            return c
    return None


@st.cache_data(show_spinner="Structure AlphaFold…")
def _fetch_af_pdb(uniprot):
    """Télécharge la structure AlphaFold d'un UniProt (b-factor = pLDDT).
    Passe par l'API pour obtenir l'URL exacte (version courante v6…)."""
    import json
    try:
        with urllib.request.urlopen(
                f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot}",
                timeout=25) as r:
            meta = json.load(r)
        url = meta[0].get("pdbUrl")
        if url:
            with urllib.request.urlopen(url, timeout=30) as r:
                return r.read().decode()
    except Exception:
        pass
    return None


def _set_bfactor(pdb_text, resmap, default=0.0):
    """Réécrit la colonne b-factor d'un PDB avec resmap {resnum: valeur}."""
    out = []
    for ln in pdb_text.splitlines():
        if ln.startswith(("ATOM", "HETATM")) and len(ln) >= 66:
            try:
                resi = int(ln[22:26])
            except ValueError:
                out.append(ln)
                continue
            out.append(ln[:60] + f"{resmap.get(resi, default):6.2f}" + ln[66:])
        else:
            out.append(ln)
    return "\n".join(out)


def _bfactor_range(pdb_text):
    vals = []
    for ln in pdb_text.splitlines():
        if ln.startswith(("ATOM", "HETATM")) and len(ln) >= 66:
            try:
                vals.append(float(ln[60:66]))
            except ValueError:
                pass
    return (min(vals), max(vals)) if vals else (0.0, 1.0)


@st.cache_data(show_spinner=False)
def structure_for(slug, uniprot, _mtime):
    """Structure 3D de l'ABP + b-factor à colorer.
    Priorité : PDB du dossier ProteoCast (sensibilité) > AlphaFold + mapping du
    Variant_score si CSV présent > AlphaFold seul (pLDDT). Renvoie dict ou None."""
    # 1) PDB fourni par ProteoCast (b-factor déjà = métrique)
    d = _result_dir_or_csv(slug)
    if d:
        pdbs = [f for f in os.listdir(d) if f.lower().endswith(".pdb")]
        pref = ([f for f in pdbs if "sensitiv" in f.lower()]
                or [f for f in pdbs if "resclass" in f.lower()] or pdbs)
        if pref:
            txt = open(os.path.join(d, pref[0])).read()
            vmin, vmax = _bfactor_range(txt)
            return {"pdb": txt, "by": "ProteoCast sensitivity (b-factor of the result)",
                    "vmin": vmin, "vmax": vmax, "metric": True}
    # 2) AlphaFold + mapping du Variant_score si CSV présent
    if not uniprot:
        return None
    af = _fetch_af_pdb(uniprot)
    if af is None:
        return None
    csv = _proteocast_csv(slug)
    if csv:
        df = pd.read_csv(csv)
        if {"Residue", "Variant_score"} <= set(df.columns):
            # sensibilité par position = -moyenne (positif, plus haut = plus sensible)
            g = df.groupby("Residue")["Variant_score"].mean()
            resmap = {int(k): -float(v) for k, v in g.items()}
            txt = _set_bfactor(af, resmap)
            vmin, vmax = _bfactor_range(txt)
            return {"pdb": txt, "by": "ProteoCast mutational sensitivity "
                    "(redder = more sensitive)", "vmin": vmin, "vmax": vmax,
                    "metric": True}
    # 3) AlphaFold seul (pLDDT) — utile en attendant ProteoCast
    vmin, vmax = _bfactor_range(af)
    return {"pdb": af, "by": "AlphaFold pLDDT confidence (ProteoCast not provided yet)",
            "vmin": vmin, "vmax": vmax, "metric": False}
