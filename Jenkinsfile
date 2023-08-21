def projectName = "solid_dmft" /* set to app/repo name */

def dockerName = projectName.toLowerCase();
/* which platform to build documentation on */
def documentationPlatform = "ubuntu-clang"
/* depend on triqs upstream branch/project */
def triqsBranch = env.CHANGE_TARGET ?: env.BRANCH_NAME
def triqsProject = '/TRIQS/triqs/' + triqsBranch.replaceAll('/', '%2F')
def cthybProject = '/TRIQS/cthyb/' + triqsBranch.replaceAll('/', '%2F')
/* whether to keep and publish the results */
def keepInstall = !env.BRANCH_NAME.startsWith("PR-")

properties([
  disableConcurrentBuilds(),
  buildDiscarder(logRotator(numToKeepStr: '10', daysToKeepStr: '30')),
  pipelineTriggers(keepInstall ? [
    upstream(
      threshold: 'SUCCESS',
      upstreamProjects: cthybProject
    )
  ] : [])
])

/* map of all builds to run, populated below */
def platforms = [:]

/****************** linux builds (in docker) */
/* Each platform must have a cooresponding Dockerfile.PLATFORM in triqs/packaging */
def dockerPlatforms = ["ubuntu-clang", "ubuntu-gcc"]
/* .each is currently broken in jenkins */
for (int i = 0; i < dockerPlatforms.size(); i++) {
  def platform = dockerPlatforms[i]
  platforms[platform] = { -> node('linux && docker && triqs') {
    stage(platform) { timeout(time: 1, unit: 'HOURS') { ansiColor('xterm') {
      checkout scm
      /* construct a Dockerfile for this base */
      sh """
      ( echo "FROM flatironinstitute/triqs:${triqsBranch}-${env.STAGE_NAME}" ; sed '0,/^FROM /d' Docker/jenkins_ci_dockerfile ) > Dockerfile.jenkins
        mv -f Dockerfile.jenkins Dockerfile
      """
      /* build and tag */
      def args = ''
      if (platform == documentationPlatform)
        args = '-DBuild_Documentation=0'
      else if (platform == "sanitize")
        args = '-DASAN=ON -DUBSAN=ON'
      def img = docker.build("flatironinstitute/${dockerName}:${env.BRANCH_NAME}-${env.STAGE_NAME}", "--build-arg APPNAME=${projectName} --build-arg BUILD_ID=${env.BUILD_TAG} --build-arg CMAKE_ARGS='${args}' .")
      catchError(buildResult: 'UNSTABLE', stageResult: 'UNSTABLE') {
        img.inside() {
          sh "make -C \$BUILD/${projectName} test CTEST_OUTPUT_ON_FAILURE=1"
        }
      }
      if (!keepInstall) {
        sh "docker rmi --no-prune ${img.imageName()}"
      }
    } } }
  } }
}

/****************** wrap-up */
def error = null
try {
  parallel platforms
  if (false/*keepInstall*/) { node('linux && docker && triqs') {
    /* Publish results */
    stage("publish") { timeout(time: 5, unit: 'MINUTES') {
      def commit = sh(returnStdout: true, script: "git rev-parse HEAD").trim()
      def release = env.BRANCH_NAME == "master" || env.BRANCH_NAME == "unstable" || sh(returnStdout: true, script: "git describe --exact-match HEAD || true").trim()
      def workDir = pwd(tmp:true)
      lock('triqs_publish') {
      /* Update documention on gh-pages branch */
      dir("$workDir/gh-pages") {
        def subdir = "${projectName}/${env.BRANCH_NAME}"
        git(url: "ssh://git@github.com/TRIQS/TRIQS.github.io.git", branch: "master", credentialsId: "ssh", changelog: false)
        sh "rm -rf ${subdir}"
        docker.image("flatironinstitute/${dockerName}:${env.BRANCH_NAME}-${documentationPlatform}").inside() {
          sh """#!/bin/bash -ex
            base=\$INSTALL/share/doc
            dir="${projectName}"
            [[ -d \$base/triqs_\$dir ]] && dir=triqs_\$dir || [[ -d \$base/\$dir ]]
            cp -rp \$base/\$dir ${subdir}
          """
        }
        sh "git add -A ${subdir}"
        sh """
          git commit --author='Flatiron Jenkins <jenkins@flatironinstitute.org>' --allow-empty -m 'Generated documentation for ${subdir}' -m '${env.BUILD_TAG} ${commit}'
        """
        // note: credentials used above don't work (need JENKINS-28335)
        sh "git push origin master"
      }
      /* Update packaging repo submodule */
      if (release) { dir("$workDir/packaging") { try {
        git(url: "ssh://git@github.com/TRIQS/packaging.git", branch: env.BRANCH_NAME, credentialsId: "ssh", changelog: false)
        // note: credentials used above don't work (need JENKINS-28335)
        sh """#!/bin/bash -ex
          dir="${projectName}"
          [[ -d triqs_\$dir ]] && dir=triqs_\$dir || [[ -d \$dir ]]
          echo "160000 commit ${commit}\t\$dir" | git update-index --index-info
          git commit --author='Flatiron Jenkins <jenkins@flatironinstitute.org>' -m 'Autoupdate ${projectName}' -m '${env.BUILD_TAG}'
          git push origin ${env.BRANCH_NAME}
        """
      } catch (err) {
        /* Ignore, non-critical -- might not exist on this branch */
        echo "Failed to update packaging repo"
      } } }
      }
    } }
  } }
} catch (err) {
  error = err
} finally {
  /* send email on build failure (declarative pipeline's post section would work better) */
  if ((error != null || currentBuild.currentResult != 'SUCCESS') && env.BRANCH_NAME != "jenkins") emailext(
    subject: "\$PROJECT_NAME - Build # \$BUILD_NUMBER - FAILED",
    body: """\$PROJECT_NAME - Build # \$BUILD_NUMBER - FAILED

Check console output at \$BUILD_URL to view full results.

Building \$BRANCH_NAME for \$CAUSE
\$JOB_DESCRIPTION

Changes:
\$CHANGES

End of build log:
\${BUILD_LOG,maxLines=60}
    """,
    to: 'nwentzell@flatironinstitute.org, dsimon@flatironinstitute.org',
    recipientProviders: [
      [$class: 'DevelopersRecipientProvider'],
    ],
    replyTo: '$DEFAULT_REPLYTO'
  )
  if (error != null) throw error
}
