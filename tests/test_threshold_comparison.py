import sys
import unittest
from pathlib import Path
import pandas as pd
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'script'))
from threshold_comparison import compare_thresholds

class ThresholdTests(unittest.TestCase):
    def test_requires_connected_actins_and_direct_partner(self):
        rows=[]
        for a,b in [('A','B'),('B','C'),('D','E')]:
            rows.append(['p',a,'Actin',b,'Actin','homo'])
        rows += [['p','A','Actin','x','Cofilin','hetero'],
                 ['p','x','Cofilin','y','Indirect','hetero']]
        entries=pd.DataFrame(rows,columns=['pdb_id','Interactor 1','Interactor 1 title',
                                          'Interactor 2','Interactor 2 title','Interface type'])
        summary=pd.DataFrame({'Interface type':['homo'],'Expect value':[0.0],'Result protein':['Actin']})
        result,_=compare_thresholds(entries,summary)
        self.assertEqual(result[0]['pdb_count'],1)
        self.assertEqual(result[1]['pdb_count'],0)  # five total, but not connected
        self.assertEqual(result[0]['additional_partner_names'],['Cofilin'])
