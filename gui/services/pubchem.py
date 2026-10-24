"""Web lookup services to resolve compound names to SMILES using PubChem, NCBI, and Cactus."""

import re
import json
import urllib.request
import urllib.parse
from typing import Optional
import streamlit as st


@st.cache_data(show_spinner=False, ttl=3600)
def fetch_pubchem_smiles(query: str) -> Optional[str]:
    """Queries PubChem, NCBI E-utilities, and Cactus to resolve compound names (Italian or English) to canonical SMILES."""
    query_clean = query.strip()
    if not query_clean:
        return None

    candidates = [query_clean]
    low = query_clean.lower()

    # Italian common organic names dictionary
    it_dict = {
        'aspirina': 'aspirin',
        'paracetamolo': 'paracetamol',
        'tachipirina': 'paracetamol',
        'glucosio': 'glucose',
        'fruttosio': 'fructose',
        'saccarosio': 'sucrose',
        'galattosio': 'galactose',
        'lattosio': 'lactose',
        'maltosio': 'maltose',
        'acqua': 'water',
        'ammoniaca': 'ammonia',
        'acido acetilsalicilico': 'acetylsalicylic acid',
        'acido ascorbico': 'ascorbic acid',
        'acido citrico': 'citric acid',
        'acido ossalico': 'oxalic acid',
        'acido tartarico': 'tartaric acid',
        'acido piruvico': 'pyruvic acid',
        'acido succinico': 'succinic acid',
        'acido malico': 'malic acid',
        'acido fumarico': 'fumaric acid',
        'acido benzoico': 'benzoic acid',
        'acido acetico': 'acetic acid',
        'acido formico': 'formic acid',
        'acido propionico': 'propionic acid',
        'acido butirrico': 'butyric acid',
        'acido lattico': 'lactic acid',
        'acido valerico': 'valeric acid',
    }
    if low in it_dict:
        candidates.append(it_dict[low])

    # Rule: 'acido ...ico' -> '...ic acid'
    m = re.match(r'^acido\s+(.+)$', low)
    if m:
        core = m.group(1)
        if core.endswith('ico'):
            candidates.append(core[:-3] + 'ic acid')
        elif core.endswith('oso'):
            candidates.append(core[:-3] + 'ous acid')

    # Italian suffixes -> English
    if low.endswith('olo'):
        candidates.append(low[:-3] + 'ol')
    elif low.endswith('ano'):
        candidates.append(low[:-3] + 'ane')
    elif low.endswith('ina'):
        candidates.append(low[:-3] + 'ine')
    elif low.endswith('osio'):
        candidates.append(low[:-4] + 'ose')

    candidates = list(dict.fromkeys(candidates))
    headers = {"User-Agent": "Energyinator-App/2.0 (Quantum Molecular Modeling Studio)"}

    # Level 1: PubChem PUG REST
    for cand in candidates:
        try:
            q = urllib.parse.quote(cand)
            url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{q}/property/CanonicalSMILES,ConnectivitySMILES,IsomericSMILES/JSON"
            req = urllib.request.Request(url, headers=headers)
            with urllib.request.urlopen(req, timeout=4) as response:
                data = json.loads(response.read().decode("utf-8"))
                props = data["PropertyTable"]["Properties"][0]
                smi = props.get("IsomericSMILES") or props.get("SMILES") or props.get("CanonicalSMILES") or props.get("ConnectivitySMILES")
                if smi:
                    return smi
        except Exception:
            continue

    # Level 2: NCBI E-Utilities
    for cand in candidates[:2]:
        try:
            q = urllib.parse.quote(cand)
            url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pccompound&term={q}&retmode=json"
            req = urllib.request.Request(url, headers=headers)
            with urllib.request.urlopen(req, timeout=4) as response:
                data = json.loads(response.read().decode("utf-8"))
                id_list = data.get("esearchresult", {}).get("idlist", [])
                if id_list:
                    cid = id_list[0]
                    u2 = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/CanonicalSMILES,ConnectivitySMILES,IsomericSMILES/JSON"
                    r2 = urllib.request.Request(u2, headers=headers)
                    with urllib.request.urlopen(r2, timeout=4) as res2:
                        d2 = json.loads(res2.read().decode("utf-8"))
                        props = d2["PropertyTable"]["Properties"][0]
                        smi = props.get("IsomericSMILES") or props.get("SMILES") or props.get("CanonicalSMILES") or props.get("ConnectivitySMILES")
                        if smi:
                            return smi
        except Exception:
            continue

    # Level 3: NIH Cactus Resolver
    for cand in candidates[:2]:
        try:
            q = urllib.parse.quote(cand)
            url = f"https://cactus.nci.nih.gov/chemical/structure/{q}/smiles"
            req = urllib.request.Request(url, headers=headers)
            with urllib.request.urlopen(req, timeout=4) as response:
                smi = response.read().decode("utf-8").strip()
                if smi and "<html" not in smi.lower() and "error" not in smi.lower():
                    return smi
        except Exception:
            continue

    return None

