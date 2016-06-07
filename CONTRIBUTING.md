# Contributing to Nek5000 

NekBox on Github follows a forking workflow.
If you are unfamiliar with forking, take a look at [this great guide](https://guides.github.com/activities/forking/).
For a more detailed description of forking workflow, here's a [slightly longer read](https://www.atlassian.com/git/tutorials/comparing-workflows/forking-workflow).

Please branch off of and open pull requests to the `develop` branch.
The `master` branch is reserved for releases.

## Porting to the new layout

If you had previously forked the Nek5000/NekBox repository and made changes, you need to match the new layout before pulling.
Fortunately, this is easy!

```
cd <PATH_TO_NEKBOX>
git filter-branch --force --index-filter \
  'git ls-files -s | sed "s-\t-&core/-" |
  GIT_INDEX_FILE=$GIT_INDEX_FILE.new \
  git update-index --index-info &&
  mv $GIT_INDEX_FILE.new $GIT_INDEX_FILE' HEAD

git filter-branch --force --tree-filter \
  'test -d core/jl && mv core/jl . || echo "Nothing"' HEAD
```

Now you should be able to pull in the rest of the new layout via pull request or the command line.
