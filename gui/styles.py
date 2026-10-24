"""CSS Styles and layout theme configuration for Energyinator GUI."""

import streamlit as st


def apply_theme():
    """Applies modern scientific design styling to the Streamlit app."""
    st.markdown("""
    <style>
        .stApp {
            font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif;
        }
        .main-title {
            font-size: 2.1rem;
            font-weight: 800;
            letter-spacing: -0.03em;
            background: linear-gradient(135deg, #3b82f6 0%, #6366f1 50%, #8b5cf6 100%);
            -webkit-background-clip: text;
            -webkit-text-fill-color: transparent;
            margin-bottom: 2px;
        }
        .sub-title {
            font-size: 0.96rem;
            color: #64748b;
            font-weight: 400;
            margin-bottom: 16px;
        }
        .info-card {
            background: rgba(255, 255, 255, 0.03);
            border: 1px solid rgba(226, 232, 240, 0.25);
            border-radius: 10px;
            padding: 14px 18px;
            margin-bottom: 14px;
            box-shadow: 0 1px 3px rgba(0, 0, 0, 0.05);
        }
        .metric-badge {
            display: inline-block;
            padding: 5px 12px;
            border-radius: 6px;
            font-size: 0.82rem;
            font-weight: 600;
            background: rgba(59, 130, 246, 0.12);
            color: #2563eb;
            border: 1px solid rgba(59, 130, 246, 0.28);
            margin-right: 6px;
            margin-bottom: 6px;
        }
        .hero-summary-card {
            background: linear-gradient(135deg, rgba(30, 41, 59, 0.95) 0%, rgba(15, 23, 42, 0.98) 100%);
            border: 1px solid rgba(59, 130, 246, 0.35);
            border-radius: 12px;
            padding: 20px 24px;
            margin: 16px 0 24px 0;
            color: #f8fafc;
            box-shadow: 0 10px 25px -5px rgba(15, 23, 42, 0.25), 0 8px 10px -6px rgba(15, 23, 42, 0.25);
        }
        .hero-title {
            font-size: 1.15rem;
            font-weight: 700;
            letter-spacing: -0.01em;
            color: #ffffff;
            margin-bottom: 14px;
            display: flex;
            align-items: center;
            justify-content: space-between;
        }
        .kpi-grid {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(160px, 1fr));
            gap: 14px;
        }
        .kpi-tile {
            background: rgba(255, 255, 255, 0.04);
            border: 1px solid rgba(255, 255, 255, 0.08);
            border-radius: 8px;
            padding: 12px 14px;
        }
        .kpi-tile-highlight {
            background: rgba(16, 185, 129, 0.12);
            border: 1px solid rgba(16, 185, 129, 0.35);
            border-radius: 8px;
            padding: 12px 14px;
        }
        .kpi-label {
            font-size: 0.75rem;
            text-transform: uppercase;
            letter-spacing: 0.05em;
            color: #94a3b8;
            font-weight: 600;
            margin-bottom: 4px;
        }
        .kpi-val {
            font-size: 1.25rem;
            font-weight: 700;
            color: #ffffff;
            line-height: 1.2;
        }
        .kpi-val-green {
            font-size: 1.25rem;
            font-weight: 700;
            color: #34d399;
            line-height: 1.2;
        }
        .kpi-sub {
            font-size: 0.78rem;
            color: #cbd5e1;
            margin-top: 2px;
        }
        .chip-container {
            display: flex;
            flex-wrap: wrap;
            gap: 8px;
            margin: 8px 0 14px 0;
        }
        .stTabs [data-baseweb="tab-list"] {
            gap: 6px;
            border-bottom: 1px solid rgba(226, 232, 240, 0.4);
            margin-bottom: 18px;
        }
        .stTabs [data-baseweb="tab"] {
            padding: 10px 18px;
            font-weight: 600;
            font-size: 0.92rem;
            border-radius: 8px 8px 0 0;
            transition: all 0.2s ease;
        }
        .stButton>button {
            border-radius: 8px;
            font-weight: 600;
            font-size: 0.9rem;
            letter-spacing: 0.01em;
            transition: all 0.2s ease;
        }
    </style>
    """, unsafe_allow_html=True)


def st_image_stretched(image, caption=None):
    """Renders images full width, compatible across Streamlit versions."""
    try:
        st.image(image, caption=caption, width="stretch")
    except (TypeError, ValueError):
        st.image(image, caption=caption, use_container_width=True)


def st_plotly_stretched(fig):
    """Renders Plotly charts full width, compatible across Streamlit versions."""
    try:
        st.plotly_chart(fig, width="stretch")
    except (TypeError, ValueError):
        st.plotly_chart(fig, use_container_width=True)


def st_dataframe_stretched(df, hide_index=True):
    """Renders DataFrames full width, compatible across Streamlit versions."""
    try:
        st.dataframe(df, width="stretch", hide_index=hide_index)
    except (TypeError, ValueError):
        st.dataframe(df, use_container_width=True, hide_index=hide_index)
