
```
git init
git add .
git commit -m "Initial commit"
git branch -M main
git remote add origin https://github.com/YOUR-USERNAME/YOUR-REPO-NAME.git
git push -u origin main
```

For subsequent commits
```
git add .
git commit -m "Add DEG overlap analysis for fibroblasts"
git push
```

For new branches
```
git checkout -b pdf-figures
git add .
git commit -m "message"
git push --set-upstream origin pdf-figures
```