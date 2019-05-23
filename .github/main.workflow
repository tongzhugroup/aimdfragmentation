workflow "Release to pypi" {
  on = "release"
  resolves = ["publish"]
}

action "publish" {
  uses = "mariamrf/py-package-publish-action@master"
  args = "check"
  secrets = ["TWINE_PASSWORD", "TWINE_USERNAME"]
  env = {
    PYTHON_VERSION = "3.7.3"
  }
}

workflow "Test" {
  on = "push"
  resolves = ["Tox"]
}

action "Tox" {
  uses = "njzjz/actions/tox-conda@master"
  secrets = [
    "COVERALLS_REPO_TOKEN",
    "CODECOV_TOKEN",
    "GAUSSIANURL",
  ]
}
