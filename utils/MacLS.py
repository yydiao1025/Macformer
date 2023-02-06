# -*- coding: utf-8 -*-

import rdkit
from rdkit import Chem
import numpy as np
#from rdkit.Geometry.rdGeometry import ComputeDihedralAngle
from rdkit.Chem.rdMolTransforms import GetDihedralDeg
from RMSD import centroid,kabsch_rmsd,kabsch,reorder_hungarian,check_reflections
import copy

def compute_distance_angle_frag(SDfile,outfile):
    frags=Chem.SDMolSupplier(SDfile)
    fw=Chem.SDWriter(outfile)

    for frag in frags:
        pos= frag.GetConformer()
        dummy_atoms=[atom for atom in frag.GetAtoms() if atom.GetSymbol()=='*']
        dummy_index=[atom.GetIdx() for atom in dummy_atoms]
        nei_index=[atom.GetNeighbors()[0].GetIdx() for atom in dummy_atoms]
        
        dummy_coords=[pos.GetAtomPosition(idx) for idx in dummy_index]
        nei_coords=[pos.GetAtomPosition(idx) for idx in nei_index]
        
        dummy_distance= np.linalg.norm(dummy_coords[0]-dummy_coords[1])
        nei_distance = np.linalg.norm(nei_coords[0]-nei_coords[1])
        
        angle=GetDihedralDeg(pos, dummy_index[0],nei_index[0],nei_index[1],dummy_index[1])
        
        frag.SetProp("_dummy_distance",str(round(dummy_distance,2)))
        frag.SetProp("_nei_distance",str(round(nei_distance,2)))
        frag.SetProp("_angle",str(round(angle,2)))
        

        fw.write(frag)
    fw.close()

def maclink(frags,linkers,outfile):
    for frag in frags:
#        n+=1
#        if n%50==0:
#            print (n)
        macros=[]
        rmsd=[]
        dummy_distance_frag=float(frag.GetProp("_dummy_distance"))
        nei_distance_frag=float(frag.GetProp("_nei_distance"))
        angle_frag=float(frag.GetProp("_angle"))
    
        dummy_frag=[a for a in frag.GetAtoms() if a.GetSymbol()=='*']

        p_view=[dummy_frag[0].GetIdx(),dummy_frag[0].GetNeighbors()[0].GetIdx(),
                dummy_frag[1].GetNeighbors()[0].GetIdx(),dummy_frag[1].GetIdx()]
    
        p_all=frag.GetConformer().GetPositions()
    
        p_coord_original=p_all[p_view]
        

    
        for linker in linkers:
            dummy_distance_linker=float(linker.GetProp("_dummy_distance"))
            nei_distance_linker=float(linker.GetProp("_nei_distance"))
            angle_linker=float(linker.GetProp("_angle"))
        
            if abs(dummy_distance_frag-nei_distance_linker)>0.5: continue
            if abs(nei_distance_frag-dummy_distance_linker)>0.5: continue
            if abs(angle_frag-angle_linker)>20: continue
        
            dummy_linker=[a for a in linker.GetAtoms() if a.GetSymbol()=='*']
#            dummy_linker=[dummy_linker2[1],dummy_linker2[0]]
        
            q_all=linker.GetConformer().GetPositions()

            q_view=[dummy_linker[0].GetNeighbors()[0].GetIdx(),dummy_linker[0].GetIdx(),
                    dummy_linker[1].GetIdx(),dummy_linker[1].GetNeighbors()[0].GetIdx()]
    
            q_atoms=np.array([linker.GetAtomWithIdx(q_view[0]).GetSymbol(),
                             frag.GetAtomWithIdx(p_view[1]).GetSymbol(),
                             frag.GetAtomWithIdx(p_view[2]).GetSymbol(),
                             linker.GetAtomWithIdx(q_view[-1]).GetSymbol()])
        
            q_coord = copy.deepcopy(q_all[q_view])
            q_cent = centroid(q_coord)
        
            p_atoms=copy.deepcopy(q_atoms)
            p_coord = copy.deepcopy(p_coord_original)
            p_cent = centroid(p_coord)
            p_coord -= p_cent
        
            result_rmsd, _, _, q_review = check_reflections(
                     p_atoms,
                     q_atoms,
                     p_coord,
                     q_coord,
                     reorder_method=None,
                     rotation_method=kabsch_rmsd,
                     keep_stereo=True,
                     )
               
            U = kabsch(q_coord, p_coord)

                # recenter all atoms and rotate all atoms
            q_all -= q_cent
            q_all = np.dot(q_all, U)

                # center q on p's original coordinates
            q_all += p_cent
        
            if result_rmsd>3.0: continue
    
            for i in range(linker.GetNumAtoms()):
                linker.GetConformer().SetAtomPosition(i,q_all[i]) 
                
            temp_lig=Chem.RWMol(frag)
            for i in range(linker.GetNumAtoms()):
                j=temp_lig.AddAtom(linker.GetAtomWithIdx(i))
                temp_lig.GetConformer().SetAtomPosition(j,q_all[i])
                
            for bs in linker.GetBonds():
                x,y=bs.GetBeginAtomIdx(),bs.GetEndAtomIdx()
                temp_lig.AddBond(x+frag.GetNumAtoms(),y+frag.GetNumAtoms(),bs.GetBondType())
            
            temp_lig.AddBond(p_view[1],q_view[0]+frag.GetNumAtoms(),Chem.BondType.SINGLE)
            temp_lig.AddBond(p_view[2],q_view[3]+frag.GetNumAtoms(),Chem.BondType.SINGLE)
            
            remove=[a.GetIdx() for a in temp_lig.GetAtoms() if a.GetSymbol()=='*']
            remove.sort(reverse=True)
            for a in remove:
                temp_lig.RemoveAtom(a) 
            temp_lig.SetProp("RMSD",str(result_rmsd))
                    
            macros.append(copy.deepcopy(temp_lig))
            rmsd.append(result_rmsd)
            
        if len(macros)==0:continue    

        macros_rmsd=zip(macros,rmsd)
        sorted_macros_rmsd=sorted(macros_rmsd,key=lambda x:x[1])
    
        sorted_macros_rmsd_split=zip(*sorted_macros_rmsd)
#    if ValueError: continue
        sorted_macros,sort_rmsd=[list(x) for x in sorted_macros_rmsd_split]
        for m in sorted_macros[:10]:
            outfile.write(m)
                  
    outfile.close()
