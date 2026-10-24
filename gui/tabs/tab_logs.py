"""Tab 5: Calculation and Quantum Engine Execution Logs."""

import os
import streamlit as st
from energyinator.models import CalculationResult


def render_tab_logs(res: CalculationResult):
    """Renders the execution run.log and raw Psi4 output."""
    st.markdown("#### Log di Esecuzione")
    log_file = os.path.join(res.output_dir, "run.log")
    if os.path.exists(log_file):
        with open(log_file, "r", encoding="utf-8") as f:
            st.code(f.read(), language="log")
    else:
        st.info("Nessun file di log trovato.")

    psi4_out = os.path.join(res.output_dir, "psi4_output.dat")
    if os.path.exists(psi4_out):
        with st.expander("Output scientifico grezzo di Psi4 (ultime righe)"):
            with open(psi4_out, "r", encoding="utf-8") as f:
                st.code(f.read()[-6000:], language="text")

