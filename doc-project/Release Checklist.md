# Neuropixel Utils Release Checklist

* Update and double-check the version number and date in `VERSION`
* Commit all the files changed by the above steps.
  * Use form: `git commit -a -m "Cut release v<version>"`
* Create a git tag and push it and the changes to GitHub.
  * `git tag v<version>`
  * `git push; git push --tags`
* Create a new GitHub release from the tag.
  * Just use `<version>` as the name for the release.
  * Upload the dist tarball as a file for the release.
* Open development for next version
  * Update version number in `VERSION` to next patch or minor version, as appropriate, and put a `-SNAPSHOT` suffix at the end of it.

* If there were any problems following these instructions exactly as written, report it as a bug.
