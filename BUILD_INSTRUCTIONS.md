# Building the PlanetProfile Python Package

This document provides complete instructions for building the PlanetProfile Python package for distribution.

## Prerequisites

Before building the package, ensure you have the following installed:

1. **Python 3.8 or higher** (Python 3.11 recommended for developers)
2. **pip** (latest version recommended)
3. **Required build tools**:
   ```bash
   pip install --upgrade pip build twine
   ```

## Package Structure

The package is configured using modern Python packaging standards with:

- **`pyproject.toml`**: Modern build configuration file (PEP 517/518 compliant)
- **`setup.py`**: Legacy setup file (maintained for backwards compatibility)
- **`MANIFEST.in`**: Specifies additional files to include in the distribution

## Building the Package

### Step 1: Clean Previous Builds (Optional but Recommended)

Before building a new version, clean up any previous build artifacts:

```bash
rm -rf dist/ build/ *.egg-info/
```

### Step 2: Build the Package

From the repository root directory, run:

```bash
python -m build
```

This command will:
- Create a `dist/` directory
- Generate two distribution files:
  - **Source distribution** (`.tar.gz`): Contains source code
  - **Wheel distribution** (`.whl`): Pre-built binary distribution

**Expected output:**
```
dist/
├── planetprofile-3.0.0-py3-none-any.whl
└── planetprofile-3.0.0.tar.gz
```

### Step 3: Verify the Build

Check that the package files were created successfully:

```bash
ls -lh dist/
```

You should see both the `.whl` and `.tar.gz` files.

## Testing the Built Package

### Local Installation Test

Before uploading to PyPI, test the package locally:

```bash
# Create a test virtual environment
python -m venv test_env
source test_env/bin/activate  # On Windows: test_env\Scripts\activate

# Install from the local wheel file
pip install dist/planetprofile-3.0.0-py3-none-any.whl

# Test the installation
python -c "import PlanetProfile; print(PlanetProfile.__file__)"

# Deactivate and remove test environment
deactivate
rm -rf test_env
```

### Check Package Contents

Inspect the package contents without installing:

```bash
# For wheel file
python -m zipfile -l dist/planetprofile-3.0.0-py3-none-any.whl

# For source distribution
tar -tzf dist/planetprofile-3.0.0.tar.gz | head -20
```

## Uploading to PyPI

### Prerequisites for Upload

1. **PyPI Account**: Create an account at https://pypi.org/account/register/
2. **API Token**: Generate an API token from your PyPI account settings
3. **Twine Installed**: Already installed in the prerequisites step

### Upload Process

```bash
# Upload to PyPI
python -m twine upload dist/* --verbose
```

When prompted:
- **Username**: Enter `__token__` (exactly as shown)
- **Password**: Paste your PyPI API token

### Upload to Test PyPI (Recommended First)

Before uploading to the real PyPI, test on Test PyPI:

```bash
# Upload to Test PyPI
python -m twine upload --repository testpypi dist/* --verbose
```

Then test installing from Test PyPI:

```bash
pip install --index-url https://test.pypi.org/simple/ PlanetProfile
```

## Updating the Version

Before building a new release:

1. **Update version number** in both:
   - `pyproject.toml` (line 6: `version = "X.Y.Z"`)
   - `PlanetProfile/Utilities/PPverNum.txt`
   
2. **Update CHANGELOG.md** with release notes

3. **Clean old build artifacts**:
   ```bash
   rm -rf dist/ build/ *.egg-info/
   ```

4. **Build and upload** following the steps above

## Troubleshooting

### Issue: `python -m build` fails with "No module named build"

**Solution:**
```bash
pip install --upgrade build
```

### Issue: Build includes unwanted files

**Solution:**
- Check `.gitignore` to exclude build artifacts
- Update `MANIFEST.in` to control included files
- Use `global-exclude` patterns in `MANIFEST.in`

### Issue: Module not found after installation

**Solution:**
- Verify `pyproject.toml` `[tool.setuptools.packages.find]` settings
- Ensure `__init__.py` files exist in all package directories
- Check the wheel contents with: `python -m zipfile -l dist/*.whl`

### Issue: Upload fails with authentication error

**Solution:**
- Ensure username is exactly `__token__` (with double underscores)
- Verify API token has correct permissions
- Check that token hasn't expired

## Build Configuration Files

### pyproject.toml

Modern Python packaging configuration file. Key sections:
- `[build-system]`: Specifies build backend (setuptools)
- `[project]`: Package metadata (name, version, dependencies)
- `[tool.setuptools]`: Setuptools-specific configuration

### setup.py

Legacy setup file maintained for backwards compatibility. The build system uses `pyproject.toml` as primary configuration.

### MANIFEST.in

Controls which non-Python files are included:
- Python files (`.py`) are included automatically
- Data files must be explicitly listed
- Uses patterns like `recursive-include`, `include`, `global-exclude`

## Continuous Integration Notes

For automated builds in CI/CD:

```bash
# Install build dependencies
pip install --upgrade pip build twine

# Build
python -m build

# Optional: Check distribution
twine check dist/*

# Upload (using environment variables for credentials)
python -m twine upload dist/* --non-interactive
```

## Additional Resources

- Python Packaging Guide: https://packaging.python.org/
- PyPI: https://pypi.org/
- setuptools documentation: https://setuptools.pypa.io/
- PEP 517: https://www.python.org/dev/peps/pep-0517/
- PEP 518: https://www.python.org/dev/peps/pep-0518/

## Support

For questions or issues related to package building, contact:
- Dr. Steven D. Vance - steven.d.vance@jpl.nasa.gov
- Open an issue: https://github.com/vancesteven/PlanetProfile/issues
