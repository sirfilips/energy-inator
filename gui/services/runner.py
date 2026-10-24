"""Calculation process runner managing CLI subprocess execution, real-time logging, and cancellation."""

import os
import re
import sys
import subprocess
from datetime import datetime
from typing import Dict, Any
import streamlit as st
from energyinator.models import CalculationResult


def run_calculation_job(cfg: Dict[str, Any]):
    """Launches the CLI calculation in an isolated subprocess with live progress parsing and abort capability."""
    smiles = cfg["smiles"]
    if not smiles:
        st.error("Inserire una stringa SMILES valida prima di avviare il calcolo.")
        return

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    script_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    run_output_dir = os.path.join(script_dir, "output", f"result_{timestamp}")
    os.makedirs(run_output_dir, exist_ok=True)

    cmd = [
        sys.executable, "-m", "energyinator.cli",
        smiles,
        "-n", str(cfg["n_confs"]),
        "-r", str(cfg["rmsd_thresh"]),
        "--post-dft-rmsd", str(cfg["post_rmsd"]),
        "--ff", cfg["ff_choice"],
        "-w", str(cfg["energy_window"] if cfg["enable_window"] else 0.0),
        "-m", cfg["method"],
        "-b", cfg["basis"],
        "-t", str(cfg["threads"]),
        "--memory", cfg["memory"],
        "--temp", str(cfg["temp_k"]),
        "-o", run_output_dir
    ]
    if cfg["limit_top_n"] and cfg["top_n_val"]:
        cmd.extend(["--top-n", str(int(cfg["top_n_val"]))])
    if cfg.get("use_xtb"):
        cmd.append("--use-xtb")
        if cfg.get("xtb_method"):
            cmd.extend(["--xtb-method", str(cfg["xtb_method"])])
        if not cfg.get("xtb_opt", True):
            cmd.append("--xtb-sp")
        if cfg.get("xtb_window"):
            cmd.extend(["--xtb-window", str(cfg["xtb_window"])])
        if cfg.get("xtb_top_n"):
            cmd.extend(["--xtb-top-n", str(int(cfg["xtb_top_n"]))])
    if cfg["sp_only"]:
        cmd.append("--sp-only")
    if cfg["solvent_val"]:
        cmd.extend(["--solvent", cfg["solvent_val"]])
    if cfg["calc_freq"]:
        cmd.append("--frequencies")
    if cfg["override_charge"]:
        cmd.extend(["--charge", str(int(cfg["manual_charge"]))])
    if cfg["override_mult"]:
        cmd.extend(["--multiplicity", str(int(cfg["manual_mult"]))])

    prog_bar = st.progress(5, text="Inizializzazione del processo quantistico isolato...")
    log_expander = st.expander("Console di Esecuzione (Output in tempo reale)", expanded=True)
    with log_expander:
        log_box = st.empty()

    env = os.environ.copy()
    env["PYTHONUNBUFFERED"] = "1"

    proc = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
        universal_newlines=True,
        env=env
    )
    st.session_state["active_proc"] = proc

    all_lines = []
    curr_pct = 5
    curr_text = "Inizializzazione..."

    try:
        for line in iter(proc.stdout.readline, ''):
            all_lines.append(line)
            lower_line = line.lower()

            if "generazione di" in lower_line or "generati" in lower_line:
                curr_pct = 15
                curr_text = "Generazione conformeri 3D con ETKDGv3..."
                prog_bar.progress(curr_pct, text=curr_text)

            elif "ottimizzazione conformeri" in lower_line or "mmff" in lower_line or "uff" in lower_line:
                curr_pct = 25
                curr_text = f"Pre-ottimizzazione con campo di forze ({cfg['ff_choice'].upper()})..."
                prog_bar.progress(curr_pct, text=curr_text)

            elif "conformeri unici selezionati" in lower_line or "finestra energetica" in lower_line:
                curr_pct = 35
                curr_text = "Filtraggio per simmetria RMSD ed energy window..."
                prog_bar.progress(curr_pct, text=curr_text)

            elif "progresso dft:" in lower_line or "psi4 quantum" in lower_line:
                match = re.search(r'\[(\d+)/(\d+)\]', line)
                if match:
                    step = int(match.group(1))
                    total = int(match.group(2))
                    pct = 35 + int((step / total) * 45)
                    curr_pct = min(85, max(curr_pct, pct))
                    curr_text = f"Calcolo quantistico Psi4 ({cfg['method'].upper()}/{cfg['basis']}): conformero {step} di {total}..."
                    prog_bar.progress(curr_pct, text=curr_text)

            elif "deduplicazione post-dft" in lower_line or "conformeri collassati" in lower_line:
                curr_pct = 88
                curr_text = "Deduplicazione geometrica post-DFT dei minimi collassati..."
                prog_bar.progress(curr_pct, text=curr_text)

            elif "calcolo frequenze" in lower_line or "frequenze vibrazionali" in lower_line:
                curr_pct = 92
                curr_text = "Calcolo frequenze vibrazionali e correzioni termiche (Gibbs)..."
                prog_bar.progress(curr_pct, text=curr_text)

            elif "report excel" in lower_line or "esportato con successo" in lower_line:
                curr_pct = 96
                curr_text = "Analisi termodinamica di Boltzmann ed esportazione dati..."
                prog_bar.progress(curr_pct, text=curr_text)

            log_box.code("".join(all_lines[-20:]), language="log")

        proc.stdout.close()
        ret_code = proc.wait()

        if ret_code == 0:
            prog_bar.progress(100, text="Calcolo completato con successo!")
            json_file = os.path.join(run_output_dir, "result.json")
            if os.path.exists(json_file):
                res = CalculationResult.from_json(json_file)
                st.session_state.calc_result = res
                st.success(f"Calcolo terminato con successo in {res.elapsed_time:.1f} secondi.")
                st.rerun()
            else:
                st.error("Il processo è terminato ma 'result.json' non è stato trovato.")
        else:
            prog_bar.empty()
            st.error(f"Errore durante l'esecuzione del calcolo (Codice di uscita: {ret_code}).")
            st.code("".join(all_lines), language="log")

    finally:
        st.session_state["active_proc"] = None


def stop_active_calculation():
    """Terminates any active calculation process."""
    proc = st.session_state.get("active_proc")
    if proc is not None:
        try:
            proc.terminate()
            proc.kill()
            st.warning("Calcolo interrotto dall'utente.")
        except Exception as e:
            st.error(f"Errore durante interruzione: {e}")
        st.session_state["active_proc"] = None
        st.rerun()

