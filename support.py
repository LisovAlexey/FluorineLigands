import typing as tp

import pandas as pd
import numpy as np

import requests
import re
import time
from copy import deepcopy

from dataclasses import dataclass
from pathlib import Path

from collections import OrderedDict, namedtuple
from itertools import combinations

import shutil

import functools

import subprocess


ligand_download_api = "https://files.rcsb.org/ligands/download/"


def exception_handler(func, attempts = 5, timeout = 60):
    def wrapper(*args, **kwargs):
        for i in range(attempts):
            try:
                return func(*args, **kwargs)
            except ValueError:
                print(f"Exception raised, sleep for {timeout} seconds")
                time.sleep(timeout)
                next
                
    return wrapper




@exception_handler
def get_subtructure_smiles(query, from_=0, rows=100):
    url = "https://search.rcsb.org/rcsbsearch/v1/query"
    
    fluorine_select_params = {
      "query": {
        "type": "terminal",
        "service": "chemical",
        "parameters": {
          "value": query,
          "type": "descriptor",
          "descriptor_type": "SMILES",
          "match_type": "sub-struct-graph-relaxed"
        }
      },
      "request_options": {
        "pager": {
          "start": from_,
          "rows": rows
        }
      },

      "return_type": "mol_definition"
    }
    
    resp = requests.post(url=url, json=fluorine_select_params)
    
    if resp.status_code == 200:
        print(f"Query '{query}', ({from_}, {rows})  successed!")
        return resp.json()
    else:
        print(f"Query '{query}', ({from_}, {rows})  failed!")
        raise ValueError(resp.status_code)


def get_molecules_by_query(query):
    """
    query: str, SMILES string for search substructure
    from_: int, results page
    pages: int, number of pages
    
    """
    
    print("Attempt to get full number of entries:")
    pre_query = get_subtructure_smiles(query, 1, 1)
    
    full_number_of_mols = pre_query["total_count"]
    print("Full number of entries is:", full_number_of_mols)
    
    batch_size = 100
        
    df_full = None
    
    for i in range(0, full_number_of_mols, batch_size):
        
        temp_json = get_subtructure_smiles(query, i, batch_size)
        df = pd.DataFrame(temp_json["result_set"]).drop("services", axis = 1)
        
        if df_full is None:
            df_full = df
        else:
            df_full = pd.concat([df_full, df])
        
        time.sleep(30)
        
    
    return df_full



@exception_handler
def download_ligand(ligand_name: str, folder_to_save: Path, filename = None) -> Path:
    
    if not folder_to_save.exists():
            folder_to_save.mkdir(parents = True, exist_ok = True)
    
    for name in [ligand_name + '_ideal.mol2', ligand_name + ".mol2"]:
        
        if (folder_to_save / name).exists():
            print(f"{name} ligand already was saved to: {folder_to_save / name}")
            return folder_to_save / name

        resp = requests.post(ligand_download_api + name)
        
        if resp.status_code == 200:
            if filename is not None:
                name = filename
            with (folder_to_save / name).open("w") as wf:
                wf.write(resp.content.decode("utf-8"))
                return folder_to_save / name
            
    if resp.status_code == 404:
        print(resp)
        print(f"Ligand {ligand_name} not found (mol2 format) ")
        return None
    
    if resp.status_code != 200:
        raise ValueError(resp.status_code, ligand_download_api + name)
        
    print(resp)
    
    
def download_ligands(ligands_list, folder_to_download, downloads_folder):
    path_mappings = {}
    
    for ligand in ligands_list:
        
        try:
            ligand_path = download_ligand(ligand, folder_to_download)
        except ValueError:
            path_mappings[ligand] = None
        
        if ligand_path is not None:
            path_mappings[ligand] = ligand_path.relative_to(downloads_folder)
        else:
            path_mappings[ligand] = None
        
    return path_mappings

def mol2_to_smiles(path) -> str:
    
    if path is None:
        return None
    
    proc = subprocess.Popen(["obabel", "-imol2", path, "-osmi"], 
                           stdout=subprocess.PIPE)
    
    out = proc.stdout.read().decode("utf-8")
    return out.split()[0]


@exception_handler
def get_pdb_by_smiles(query, from_=0, rows=100):
    url = "https://search.rcsb.org/rcsbsearch/v1/query"
    
    fluorine_select_params = {
      "query": {
        "type": "terminal",
        "service": "chemical",
        "parameters": {
          "value": query,
          "type": "descriptor",
          "descriptor_type": "SMILES",
          "match_type": "graph-exact"
        }
      },
      "request_options": {
        "pager": {
          "start": from_,
          "rows": rows
        }
      },

      "return_type": "entry"
    }
    
    resp = requests.post(url=url, json=fluorine_select_params)
    
    if resp.status_code == 200:
        print(f"Query '{query}', ({from_}, {rows})  successed!")
        return resp.json()
    elif resp.status_code == 204:
        print("Nothing found for ", query)
        return None
    else:
        print(f"Query '{query}', ({from_}, {rows})  failed!")
        raise ValueError(resp.status_code)

        
@exception_handler
def get_pdb_by_ligand_identifier(query, from_=0, rows=100):
    url = "https://search.rcsb.org/rcsbsearch/v1/query"
    
    fluorine_select_params = {
      "query": {
        "type": "terminal",
        "service": "text_chem",
        "parameters": {
          "operator": "exact_match",
          "value": query,
          "attribute": "rcsb_chem_comp_container_identifiers.rcsb_id"
        }
      },
      "request_options": {
        "pager": {
          "start": from_,
          "rows": rows
        }
      },

      "return_type": "entry"
    }
    
    resp = requests.post(url=url, json=fluorine_select_params)
    
    if resp.status_code == 200:
        print(f"Query '{query}', ({from_}, {rows})  successed!")
        return resp.json()
    elif resp.status_code == 204:
        print("Nothing found for ", query)
        return None
    else:
        print(f"Query '{query}', ({from_}, {rows})  failed!")
        raise ValueError(resp.status_code)
        

def get_list_of_pdb_by_smiles(query):
    batch_size = 2000
    
    test_response = get_pdb_by_smiles(query, 0, batch_size)
    
    if test_response is None:
        return []
    
    num_of_entries = test_response['total_count']
    
    if num_of_entries > 4 * batch_size:
        return num_of_entries

    
    if num_of_entries < batch_size:
        return list(pd.DataFrame(test_response['result_set'])["identifier"])
    
    else:
        
        full_list = None
        
        for i in range(0, num_of_entries, batch_size):

            temp_json = get_pdb_by_smiles(query, i, batch_size)
            pdb_ids = list(pd.DataFrame(temp_json["result_set"])["identifier"])

            if full_list is None:
                full_list = []
            else:
                full_list.extend(pdb_ids)

            time.sleep(30)
            
        return full_list

def get_list_of_pdb_by_ligand_identifier(query):
    batch_size = 2000
    
    test_response = get_pdb_by_ligand_identifier(query, 0, batch_size)
    
    if test_response is None:
        return []
    
    num_of_entries = test_response['total_count']
    
    if num_of_entries < batch_size:
        return list(pd.DataFrame(test_response['result_set'])["identifier"])
    
    else:
        
        full_list = None
        
        for i in range(0, num_of_entries, batch_size):

            temp_json = get_pdb_by_ligand_identifier(query, i, batch_size)
            pdb_ids = list(pd.DataFrame(temp_json["result_set"])["identifier"])

            if full_list is None:
                full_list = []
            else:
                full_list.extend(pdb_ids)

            time.sleep(30)
            
        return full_list

    

from atoms import Atom, Mol2Atoms, Mol2TriposMolecule    
    

def do_atom_replacement(mol2file: Path, atom, atom_replacement, save_to: Path):

    mol2 = Mol2TriposMolecule.from_mol2(mol2file)
    
    # Generate one-length atom replacements
    atom_replacements = mol2.generate_atom_replacements(atom, atom_replacement, rep_lens = [1])

    if not save_to.exists():
        save_to.mkdir(exist_ok=True, parents=True)

    for i, replaced_atoms in enumerate(atom_replacements):
        mol2.atoms = replaced_atoms
        mol2.to_mol2(save_to / (mol2file.stem + f"_{i}.mol2"))
        
        
@functools.lru_cache(maxsize=32000)
def get_activity_df_by_pdb(pdb, ligand_comp_id = None):
    print(pdb)
    
    full_df = None
    entity_id = 1
    r = requests.get(f"https://data.rcsb.org/rest/v1/core/entry/{pdb}")
    js = r.json()
    
    if 'rcsb_binding_affinity' in js:
        df = pd.DataFrame(js['rcsb_binding_affinity'])

        return df
    
@functools.lru_cache()
def get_subunits_by_pdb_and_ligand_id(pdb, ligand_id, from_=0, rows=100):
    
    url = "https://search.rcsb.org/rcsbsearch/v1/query"
    
    query = f"{pdb} {ligand_id}"
    
    fluorine_select_params = {
      "query": {
    "type": "group",
    "logical_operator": "and",
    "nodes": [
      {
        "type": "terminal",
        "service": "text",
        "parameters": {
          "attribute": "rcsb_entry_container_identifiers.entry_id",
          "operator": "exact_match",
          "value": pdb
        }
      },
      {
        "type": "terminal",
        "service": "text_chem",
        "parameters": {
          "attribute": "rcsb_chem_comp_container_identifiers.comp_id",
          "operator": "exact_match",
          "value": ligand_id
        }
      }
    ]
  },
      "return_type": "polymer_entity",
      "request_options": {
        "pager": {
          "start": 0,
          "rows": 25
        }
    }
    }
    
    resp = requests.post(url=url, json=fluorine_select_params)
    
    if resp.status_code == 200:
        print(f"Query '{query}', ({from_}, {rows})  successed!")
        id_list = list(pd.DataFrame(resp.json()["result_set"])["identifier"].unique())
        return list(map(lambda x: x.split("_")[-1], id_list))
    elif resp.status_code == 204:
        print("Nothing found for ", query)
        return None
    else:
        print(f"Query '{query}', ({from_}, {rows})  failed!")
        raise ValueError(resp.status_code)
        

def get_uniprot_by_pdb_id_entity_id(pdb, entity_id):
    """
    Get UNIPROT by PDB ID and entity id
    """
    if entity_id is None:
        raise ValueError(f"entity_id is None")
        
    query = f"{pdb} {entity_id}"
    url = f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb}/{entity_id}"
    
    resp = requests.get(url=url)
    if resp.status_code == 200:
        print(f"Query '{query}'  successed!")
        if 'uniprot_ids' in resp.json()["rcsb_polymer_entity_container_identifiers"]:
            return resp.json()["rcsb_polymer_entity_container_identifiers"]['uniprot_ids']
        else:
            return []
    elif resp.status_code == 204:
        print("Nothing found for ", query)
        return None
    else:
        print(f"Query '{query}'  failed!")
        raise ValueError(resp.status_code)
        
@functools.lru_cache(maxsize = 20000)
def get_uniprot_by_pdb_id(pdb, entity_id):
    """
    Get UNIPROT by PDB ID and entity id.
    Same as get_uniprot_by_pdb_id_entity_id, but if entity_id is None 
    do scan of all uniprots for PDB
    """
    
    if entity_id is not None:
        return get_uniprot_by_pdb_id_entity_id(pdb, entity_id)
    else:
        l = []
        none_counter = 0
        for entity_id in range(1, 10):
            try:
                if none_counter == 2:
                    break
                
                resp = get_uniprot_by_pdb_id_entity_id(pdb, entity_id)
                if isinstance(resp, list):
                    l.extend(resp)
                    
            except KeyboardInterrupt:
                return
            except BaseException:
                none_counter += 1
                
        
        l = list(set(l))
        
    return l

@functools.lru_cache(maxsize = 20000)
def pdb2pubmed(pdb):
    query = pdb
    url = f"https://data.rcsb.org/rest/v1/core/pubmed/{pdb}"
    
    resp = requests.get(url=url)
    if resp.status_code == 200:
        print(f"Query '{query}'  successed!")
        if 'pubmed_id' in resp.json()["rcsb_pubmed_container_identifiers"]:
            return resp.json()["rcsb_pubmed_container_identifiers"]['pubmed_id']
        else:
            return None
    elif resp.status_code == 204 or resp.status_code == 404:
        print("Nothing found for ", query)
        return None
    else:
        print(f"Query '{query}'  failed!")
        raise ValueError(resp.status_code)
    
        