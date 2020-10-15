#https://stackoverflow.com/questions/16524225/how-can-i-populate-the-git-commit-id-into-a-file-when-i-commit/38087913#38087913
# date formatting: https://stackoverflow.com/questions/7853332/how-to-change-git-log-date-formats
git log -n 1 --date=local --format=format:"#define GIT_COMMIT \"%h (%ad) \"" HEAD > gitcommit.h
