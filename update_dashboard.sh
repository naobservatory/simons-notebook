#! /bin/bash

date=$(date '+%Y-%m-%d')

# Exit on error
set -e

# Render dashboard and push changes
quarto render dashboard/index.qmd
git add docs/dashboard/index.html docs/search.json
git commit -m "Updated dashboard on $date"
git push

# Update remote server
ssh sb-data "cd simons-notebook && git pull"
