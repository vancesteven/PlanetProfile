"""Cached first-page PDF preview for the GUI.

Every page that shows a PDF figure used to call pdf2image.convert_from_path
on every Streamlit rerun, converting every page of every PDF visible in any
tab. This module renders only the first page, once, and caches the PNG bytes
keyed on (path, mtime, dpi) so a regenerated PDF invalidates its own entry.
"""
import os
from io import BytesIO

import streamlit as st


@st.cache_data(max_entries=64, show_spinner=False)
def _pdf_first_page_png(path: str, mtime: float, dpi: int = 150) -> bytes:
    # mtime is part of the cache key: a re-run that overwrites the PDF
    # changes mtime and forces a fresh conversion.
    from pdf2image import convert_from_path
    images = convert_from_path(path, dpi=dpi, first_page=1, last_page=1)
    buf = BytesIO()
    images[0].save(buf, format="PNG")
    return buf.getvalue()


def show_pdf_preview(path: str, caption: str | None = None,
                     key: str | None = None, dpi: int = 150) -> bool:
    """Render the first page of a PDF (cached). Falls back to a download
    button when conversion fails (e.g. poppler missing). Returns True when
    the image rendered."""
    try:
        png = _pdf_first_page_png(path, os.path.getmtime(path), dpi)
        st.image(png, width="stretch", caption=caption)
        return True
    except Exception as exc:
        st.warning(f"PDF preview unavailable ({exc}); download instead:")
        try:
            with open(path, 'rb') as fpdf:
                st.download_button("⬇️ Download PDF", fpdf,
                                   file_name=os.path.basename(path),
                                   key=key or f"dl_{path}")
        except OSError as open_exc:
            st.error(f"Cannot read PDF: {open_exc}")
        return False
