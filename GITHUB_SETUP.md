# Publish XDiffCalculator to GitHub

These steps are for the first publication. They are instructions only; no repository has been initialized, committed or pushed automatically.

## 1. Create an empty GitHub repository

Sign in to GitHub and select **New repository**. Use `XDiffCalculator` as the name, choose visibility, and create it without adding a README, .gitignore or license through the wizard. The local project already has a README and .gitignore.

## 2. Initialize locally

Open PowerShell:

```powershell
cd C:\Users\hfang\Documents\GitHub\XDiffCalculator
git init -b main
git config user.name "Haixing Fang"
git config user.email "haixingfang868@gmail.com"
git status
```

The identity settings apply only to this repository.

## 3. Review and commit

```powershell
git add .
git diff --cached --stat
git status
git commit -m "Initial release of XDiffCalculator"
```

Review the staged list before committing. `.gitignore` excludes `.venv/`, `build/`, `dist/`, bytecode, generated spec files and smoke reports. Source, requirements, documentation, tests and examples should be included.

## 4. Connect and push

Use the URL shown by your new repository. For the account `haixingfang`:

```powershell
git remote add origin https://github.com/haixingfang/XDiffCalculator.git
git remote -v
git push -u origin main
```

Complete the Git authentication prompt if shown. Refresh the repository page to see the source and rendered README.

## 5. Attach the standalone program

Open **Releases → Draft a new release**, create tag `v0.1.0` targeting `main`, enter a title/description, and attach **dist/XDiffCalculator.exe**. Publish when ready. The source archive does not include the ignored executable, so users should download the release asset.

For later changes, review `git status`, add the intended files, commit and push. Rebuild before uploading a new executable under a new release version.

Official instructions: [add local code](https://docs.github.com/en/migrations/importing-source-code/using-the-command-line-to-import-source-code/adding-locally-hosted-code-to-github), [manage releases](https://docs.github.com/en/repositories/releasing-projects-on-github/managing-releases-in-a-repository).
