"""Streamlit entry point: native navigation, one scientific page at a time."""
import streamlit as st
from app_navigation import run_navigation

st.set_page_config(layout="wide", page_title="Actin–ABP analysis · PPI3D")
st.markdown(
    "<style>h2,h3{scroll-margin-top:5rem;}"
    "[data-testid='stSidebarNav'] a[aria-current='page']{font-weight:700;}"
    "[data-testid='stSidebarNavLink']{height:auto;min-height:2rem;}"
    "[data-testid='stSidebarNavLink'] span,[data-testid='stSidebarNavLink'] p"
    "{white-space:normal;overflow:visible;text-overflow:clip;}"
    ".hoverlayer path.legend3dandfriends[style*=url]{opacity:0;}"
    "</style>", unsafe_allow_html=True,
)
run_navigation()
