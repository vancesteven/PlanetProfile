"""
Session Management for PlanetProfileApp
Handles saving, loading, and managing user sessions
"""

import json
import os
import re
from datetime import datetime
from pathlib import Path
import streamlit as st


# --- Import validation (security boundary) ---------------------------------
# Session files can come from untrusted uploads on the public deployment.
# Restoring arbitrary keys into st.session_state is an injection vector:
# e.g. ChosenPlanet feeds os.path.join() on several pages (path traversal)
# and template-module resolution. Only keys matching these patterns, with
# scalar-ish values, are ever applied; everything else is dropped and
# reported. Never unpickle session files.
_ALLOWED_KEY_RE = re.compile(
    r'^(ChosenPlanet|custom_planet_name'
    r'|changed_[A-Za-z0-9_]+|custom_[A-Za-z0-9_]+'
    r'|explore_[A-Za-z0-9_ ]+|amort_[A-Za-z0-9_ .()\-]+'
    r'|inference_[A-Za-z0-9_]+|exc_[A-Za-z0-9_ ]+'
    r'|[xy]_(?:min|max|param_select|driver_select)'
    r'|z_param_select|n[xy]_input)$'
)
_MAX_STR = 500
_MAX_COLLECTION = 300


def _known_bodies():
    default_dir = Path(__file__).parent.parent.parent / 'PlanetProfile' / 'Default'
    try:
        return {d.name for d in default_dir.iterdir() if d.is_dir()}
    except OSError:
        return set()


def _scalar_ok(v):
    return v is None or isinstance(v, bool) or (
        isinstance(v, (int, float)) and abs(float(v)) < 1e30) or (
        isinstance(v, str) and len(v) <= _MAX_STR)


def validate_session_state(state):
    """Filter an untrusted session-state dict down to safe entries.

    Returns (clean_state, rejected_keys). Values may be scalars, flat lists
    of scalars, or flat str->scalar dicts. ChosenPlanet must name a known
    body. Anything else is rejected, never raised on.
    """
    clean, rejected = {}, []
    if not isinstance(state, dict):
        return {}, ['<state is not a dict>']
    for key, value in state.items():
        ok = (isinstance(key, str) and _ALLOWED_KEY_RE.match(key))
        if ok:
            if isinstance(value, list):
                ok = len(value) <= _MAX_COLLECTION and all(_scalar_ok(v) for v in value)
            elif isinstance(value, dict):
                ok = (len(value) <= _MAX_COLLECTION and
                      all(isinstance(k, str) and len(k) <= _MAX_STR and _scalar_ok(v)
                          for k, v in value.items()))
            else:
                ok = _scalar_ok(value)
        if ok and key == 'ChosenPlanet':
            ok = value in _known_bodies()
        if ok:
            clean[key] = value
        else:
            rejected.append(str(key)[:80])
    return clean, rejected


class SessionManager:
    """Manages saving and loading of user sessions"""

    def __init__(self, sessions_dir="sessions"):
        """
        Initialize session manager

        Args:
            sessions_dir: Directory to store session files (relative to app root)
        """
        # Get app directory
        app_dir = Path(__file__).parent.parent
        self.sessions_dir = app_dir / sessions_dir

        # Create sessions directory if it doesn't exist
        self.sessions_dir.mkdir(exist_ok=True)

    def save_session(self, session_name, session_state=None, metadata=None):
        """
        Save current session to file

        Args:
            session_name: Name for the session
            session_state: Streamlit session state dict (uses st.session_state if None)
            metadata: Additional metadata dict

        Returns:
            Path to saved session file
        """
        if session_state is None:
            session_state = dict(st.session_state)

        # Create session data structure
        session_data = {
            'name': session_name,
            'timestamp': datetime.now().isoformat(),
            'metadata': metadata or {},
            'state': {}
        }

        # Filter session state to only include serializable items
        for key, value in session_state.items():
            try:
                # Test if serializable
                json.dumps(value)
                session_data['state'][key] = value
            except (TypeError, ValueError):
                # Skip non-serializable items (like Planet objects)
                # Instead, save just the essential info
                if key == 'Planet' and value is not None:
                    session_data['state']['ChosenPlanet'] = session_state.get('ChosenPlanet', 'Unknown')
                elif key.startswith('changed_') or key.startswith('custom_'):
                    # Try to convert to dict
                    try:
                        session_data['state'][key] = dict(value) if hasattr(value, '__dict__') else str(value)
                    except:
                        pass

        # Save to file
        filename = self._sanitize_filename(session_name) + '.json'
        filepath = self.sessions_dir / filename

        with open(filepath, 'w') as f:
            json.dump(session_data, f, indent=2)

        return filepath

    def load_session(self, session_name_or_path):
        """
        Load session from file

        Args:
            session_name_or_path: Session name or full path to session file

        Returns:
            Session data dict
        """
        # Determine file path
        if os.path.exists(session_name_or_path):
            filepath = Path(session_name_or_path)
        else:
            filename = self._sanitize_filename(session_name_or_path) + '.json'
            filepath = self.sessions_dir / filename

        if not filepath.exists():
            raise FileNotFoundError(f"Session not found: {filepath}")

        with open(filepath, 'r') as f:
            session_data = json.load(f)

        return session_data

    def apply_session(self, session_data):
        """
        Apply loaded session data to current session state.

        The state dict is treated as UNTRUSTED regardless of source (files
        can be uploaded): only allowlisted keys with scalar-ish values are
        applied (see validate_session_state); the rest are dropped.

        Returns:
            List of rejected key names (empty when everything applied).
        """
        clean, rejected = validate_session_state(session_data.get('state', {}))
        for key, value in clean.items():
            st.session_state[key] = value

        # Set flag to reload Planet object on main settings page
        st.session_state['need_planet_reload'] = True
        return rejected

    def list_sessions(self, sort_by='timestamp', ascending=False):
        """
        List all saved sessions

        Args:
            sort_by: 'name' or 'timestamp'
            ascending: Sort order

        Returns:
            List of (name, timestamp, metadata) tuples
        """
        sessions = []

        for filepath in self.sessions_dir.glob('*.json'):
            try:
                with open(filepath, 'r') as f:
                    data = json.load(f)
                    sessions.append({
                        'name': data['name'],
                        'filename': filepath.name,
                        'timestamp': data['timestamp'],
                        'metadata': data.get('metadata', {})
                    })
            except Exception as e:
                print(f"Error loading session {filepath}: {e}")

        # Sort
        if sort_by == 'timestamp':
            sessions.sort(key=lambda x: x['timestamp'], reverse=not ascending)
        else:
            sessions.sort(key=lambda x: x['name'], reverse=not ascending)

        return sessions

    def delete_session(self, session_name):
        """
        Delete a saved session

        Args:
            session_name: Name of session to delete

        Returns:
            True if deleted, False if not found
        """
        filename = self._sanitize_filename(session_name) + '.json'
        filepath = self.sessions_dir / filename

        if filepath.exists():
            filepath.unlink()
            return True
        return False

    def get_recent_sessions(self, n=5):
        """
        Get n most recent sessions

        Args:
            n: Number of sessions to return

        Returns:
            List of recent session dicts
        """
        all_sessions = self.list_sessions(sort_by='timestamp', ascending=False)
        return all_sessions[:n]

    def export_session(self, session_name):
        """
        Export session as JSON string for download

        Args:
            session_name: Name of session to export

        Returns:
            JSON string
        """
        session_data = self.load_session(session_name)
        return json.dumps(session_data, indent=2)

    _MAX_IMPORT_BYTES = 512 * 1024

    def parse_uploaded_session(self, json_string):
        """Parse an UNTRUSTED uploaded session file into validated data.

        Size-capped JSON only (never pickle); the state dict is filtered
        through validate_session_state. Returns (session_data, rejected).
        """
        if hasattr(json_string, 'read'):
            raw = json_string.read()
        else:
            raw = json_string
        if isinstance(raw, str):
            raw = raw.encode()
        if len(raw) > self._MAX_IMPORT_BYTES:
            raise ValueError(
                f"Session file too large ({len(raw)} bytes > "
                f"{self._MAX_IMPORT_BYTES}) — not a valid configuration file.")
        session_data = json.loads(raw)
        if not isinstance(session_data, dict):
            raise ValueError("Not a session file (top level is not an object).")
        clean, rejected = validate_session_state(session_data.get('state', {}))
        session_data['state'] = clean
        session_data['name'] = str(session_data.get('name', 'imported'))[:80]
        return session_data, rejected

    def import_session(self, json_string, session_name=None):
        """
        Import session from JSON string (validated) and save server-side.

        Args:
            json_string: JSON string or file-like object (UNTRUSTED)
            session_name: Optional new name (uses original if None)

        Returns:
            Path to saved session file
        """
        session_data, _rejected = self.parse_uploaded_session(json_string)

        # Use new name if provided
        if session_name:
            session_data['name'] = session_name

        # Save it
        filename = self._sanitize_filename(session_data['name']) + '.json'
        filepath = self.sessions_dir / filename

        with open(filepath, 'w') as f:
            json.dump(session_data, f, indent=2)

        return filepath

    @staticmethod
    def _sanitize_filename(name):
        """Sanitize filename"""
        import re
        name = re.sub(r'[<>:"/\\|?*]', '_', name)
        name = name.replace(' ', '_')
        name = name.strip('_.')
        return name


# Streamlit UI components for session management

def _current_session_json(name):
    """Serialize the CURRENT session state to a downloadable JSON config
    (no server-side write — safe for shared public deployments)."""
    sm = SessionManager()
    state = {}
    for key, value in dict(st.session_state).items():
        try:
            json.dumps(value)
            state[key] = value
        except (TypeError, ValueError):
            continue
    clean, _ = validate_session_state(state)
    return json.dumps({
        'name': name,
        'timestamp': datetime.now().isoformat(),
        'metadata': {'planet': st.session_state.get('ChosenPlanet', 'Unknown')},
        'state': clean,
    }, indent=2)


def show_session_manager_ui():
    """Show session manager UI in sidebar"""
    from Utilities.app_mode import public_mode
    sm = SessionManager()

    st.sidebar.markdown("---")
    st.sidebar.subheader("💾 Session Management")

    # Configuration as a FILE (download/upload): works everywhere, and is
    # the ONLY mechanism offered publicly — the server-side sessions dir is
    # shared by all visitors of the same container (cross-visitor leakage
    # and planted-session risk), so it is hidden in public mode.
    with st.sidebar.expander("Save / load configuration file"):
        dl_name = st.text_input(
            "Configuration name",
            value=f"PPconfig_{datetime.now().strftime('%Y%m%d_%H%M')}",
            key="config_download_name")
        st.download_button(
            "📥 Download current configuration",
            data=_current_session_json(dl_name),
            file_name=f"{SessionManager._sanitize_filename(dl_name)}.json",
            mime="application/json")
        uploaded_cfg = st.file_uploader(
            "Load configuration (JSON)", type=['json'],
            key="config_upload_file",
            help="Only recognized settings are applied; anything else in "
                 "the file is ignored.")
        if uploaded_cfg is not None and st.button("📂 Apply configuration"):
            try:
                data, rejected = sm.parse_uploaded_session(uploaded_cfg)
                sm.apply_session(data)
                msg = f"✅ Applied: {data['name']}"
                if rejected:
                    msg += f" ({len(rejected)} unrecognized entries ignored)"
                st.success(msg)
                st.rerun()
            except Exception as e:
                st.error(f"❌ Invalid configuration file: {e}")

    if public_mode():
        return  # no shared server-side sessions on the public deployment

    # Save current session
    with st.sidebar.expander("Save Current Session"):
        session_name = st.text_input(
            "Session Name",
            value=f"Session_{datetime.now().strftime('%Y%m%d_%H%M')}",
            key="save_session_name"
        )

        if st.button("💾 Save Session"):
            try:
                metadata = {
                    'planet': st.session_state.get('ChosenPlanet', 'Unknown'),
                    'has_changes': any(st.session_state.get('changed_bulk_settings_flags', {}).values())
                }
                filepath = sm.save_session(session_name, metadata=metadata)
                st.success(f"✅ Session saved: {session_name}")
            except Exception as e:
                st.error(f"❌ Error saving session: {e}")

    # Load session
    with st.sidebar.expander("Load Session"):
        sessions = sm.list_sessions()

        if not sessions:
            st.info("No saved sessions")
        else:
            session_names = [s['name'] for s in sessions]
            selected = st.selectbox(
                "Select Session",
                options=session_names,
                key="load_session_select"
            )

            col1, col2 = st.columns(2)

            with col1:
                if st.button("📂 Load"):
                    try:
                        session_data = sm.load_session(selected)
                        sm.apply_session(session_data)
                        st.success(f"✅ Loaded: {selected}")
                        st.rerun()
                    except Exception as e:
                        st.error(f"❌ Error: {e}")

            with col2:
                if st.button("🗑️ Delete"):
                    if sm.delete_session(selected):
                        st.success(f"🗑️ Deleted: {selected}")
                        st.rerun()

    # Export/Import
    with st.sidebar.expander("Export/Import"):
        # Export
        sessions = sm.list_sessions()
        if sessions:
            export_name = st.selectbox(
                "Export Session",
                options=[s['name'] for s in sessions],
                key="export_session_select"
            )

            if st.button("📤 Export to JSON"):
                json_data = sm.export_session(export_name)
                st.download_button(
                    label="📥 Download JSON",
                    data=json_data,
                    file_name=f"{export_name}.json",
                    mime="application/json"
                )

        # Import
        st.markdown("**Import Session:**")
        uploaded = st.file_uploader(
            "Upload JSON",
            type=['json'],
            key="import_session_file"
        )

        if uploaded:
            try:
                filepath = sm.import_session(uploaded)
                st.success(f"✅ Imported session")
                st.rerun()
            except Exception as e:
                st.error(f"❌ Error importing: {e}")


def show_recent_sessions():
    """Show quick access to recent sessions"""
    from Utilities.app_mode import public_mode
    if public_mode():
        return  # server-side sessions are shared across visitors — hidden
    sm = SessionManager()
    recent = sm.get_recent_sessions(3)

    if recent:
        st.sidebar.markdown("---")
        st.sidebar.subheader("📋 Recent Sessions")

        for session in recent:
            timestamp = datetime.fromisoformat(session['timestamp'])
            time_str = timestamp.strftime("%m/%d %H:%M")

            if st.sidebar.button(f"📂 {session['name']}", key=f"recent_{session['filename']}"):
                session_data = sm.load_session(session['name'])
                sm.apply_session(session_data)
                st.rerun()

            st.sidebar.caption(f"   {time_str} - {session['metadata'].get('planet', 'Unknown')}")
