# Quick Start: Building PlanetProfile Package

## TL;DR - Build the Package

```bash
# Install build tools (one-time setup)
pip install --upgrade pip build twine

# Build the package
python -m build

# Result: Two files created in dist/
# - planetprofile-3.0.0-py3-none-any.whl
# - planetprofile-3.0.0.tar.gz
```

## What Was Fixed

The CONTRIBUTING.md mentioned running `python -m build`, but the necessary configuration file (`pyproject.toml`) was missing and was actually excluded in `.gitignore`. This has now been fixed.

## Key Files

1. **`pyproject.toml`** (NEW)
   - Modern Python packaging configuration
   - Contains all package metadata and build settings
   - Required for `python -m build` to work

2. **`BUILD_INSTRUCTIONS.md`** (NEW)
   - Complete step-by-step guide for building packages
   - Includes troubleshooting and testing instructions
   - Covers PyPI upload process

3. **`.gitignore`** (UPDATED)
   - Removed `pyproject.toml` from ignore list (it should be tracked)
   - Added `build/` directory to ignore list

4. **`CONTRIBUTING.md`** (UPDATED)
   - Added reference to BUILD_INSTRUCTIONS.md
   - Updated to mention `pyproject.toml` instead of just `setup.py`

## Next Steps for Building a Release

When you're ready to create a new release:

1. Update version in **both** places:
   - `pyproject.toml` (line 6)
   - `PlanetProfile/Utilities/PPverNum.txt`

2. Clean and build:
   ```bash
   rm -rf dist/ build/ *.egg-info/
   python -m build
   ```

3. Upload to PyPI:
   ```bash
   python -m twine upload dist/* --verbose
   ```
   - Username: `__token__`
   - Password: Your PyPI API token

## Testing Locally

```bash
# Create test environment
python -m venv test_env
source test_env/bin/activate

# Install from local build
pip install dist/planetprofile-3.0.0-py3-none-any.whl

# Test import
python -c "import PlanetProfile; print('Success!')"

# Clean up
deactivate
rm -rf test_env
```

## For More Details

See **BUILD_INSTRUCTIONS.md** for comprehensive documentation including:
- Prerequisites setup
- Detailed build process
- Troubleshooting guide
- PyPI upload instructions
- CI/CD integration notes
