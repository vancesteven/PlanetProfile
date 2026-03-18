"""
Session Management for PlanetProfileApp
Handles saving, loading, and managing user sessions
"""

import pickle
import json
import os
from datetime import datetime
from pathlib import Path
import streamlit as st


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
        Apply loaded session data to current session state

        Args:
            session_data: Session data dict from load_session()
        """
        # Apply state
        for key, value in session_data['state'].items():
            st.session_state[key] = value

        # Set flag to reload Planet object on main settings page
        st.session_state['need_planet_reload'] = True

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

    def import_session(self, json_string, session_name=None):
        """
        Import session from JSON string

        Args:
            json_string: JSON string or file-like object
            session_name: Optional new name (uses original if None)

        Returns:
            Path to saved session file
        """
        # Parse JSON
        if hasattr(json_string, 'read'):
            # It's a file-like object
            session_data = json.load(json_string)
        else:
            session_data = json.loads(json_string)

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

def show_session_manager_ui():
    """Show session manager UI in sidebar"""
    sm = SessionManager()

    st.sidebar.markdown("---")
    st.sidebar.subheader("💾 Session Management")

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
