# GitHub Release Instructions for QMCkl v1.1.0

## Release Title
```
QMCkl v1.1.0 - Forces Module and Enhanced Compiler Support
```

## Release Tag
```
v1.1.0
```

## Target Branch
```
master
```

## Release Description

Use the content from `RELEASE_NOTES.md` as the GitHub release description.

## Steps to Create the Release

1. **Merge this PR** to master branch
2. **Create a new tag** on the master branch:
   ```bash
   git tag -a v1.1.0 -m "Release v1.1.0"
   git push origin v1.1.0
   ```
3. **Create GitHub Release**:
   - Go to https://github.com/TREX-CoE/qmckl/releases/new
   - Tag: v1.1.0
   - Title: QMCkl v1.1.0 - Forces Module and Enhanced Compiler Support
   - Description: Copy content from RELEASE_NOTES.md
   - Upload the tarball if available

4. **Build and Upload Source Distribution** (if not done automatically):
   ```bash
   ./autogen.sh
   ./configure
   make dist
   ```
   This creates `qmckl-1.1.0.tar.gz` which should be uploaded as a release asset.

## Key Highlights to Emphasize

- **NEW Forces Module** - Major feature for molecular dynamics and geometry optimization
- **NEW Single-Electron Move Jastrow** - Performance improvement for QMC calculations
- **Modern Intel Compiler Support** - ifx and icx compilers
- **Enhanced Debugging** - Sanitizer support
- **Security** - Flawfinder integration

## Notes

- The NEWS file contains user-facing release notes
- The RELEASE_NOTES.md contains the detailed GitHub release description
- Version has been updated in configure.ac to 1.1.0
- All changes have been committed to the copilot/create-new-release branch
