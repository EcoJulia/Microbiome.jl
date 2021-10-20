# Contributing to Microbiome.jl and BiobakeryUtils.jl

:+1::tada: First off, thanks for taking the time to contribute! :tada::+1:

## Etiquette and conduct

Contributors and maintainers for Microbiome.jl and BiobakeryUtils.jl
are expected to abide by the [Julia Community Standards](https://julialang.org/community/standards/)

## How can I contribute?

[(Adapted from the BioJulia's contributing guidelines)](https://github.com/BioJulia/Contributing)
### Reporting Bugs

If you follow the advice here, we will
better understand your report :pencil:, be able to reproduce the behaviour
:computer: :computer:, and identify related problems :mag_right:.

#### Before creating a bug report:

Please do the following:

1. Check the GitHub issue list for the package that is giving you problems.

2. If you find an issue already open for your problem, add a comment to let
  everyone know that you are experiencing the same issue.

3. If no **currently open** issue already exists for your problem that has already been
   then you should create a new issue.

   > **Note:** If you find a **Closed** issue that seems like it is the same thing
   > that you're experiencing, open a new issue and include a link to the original
   > issue in the body of your new one.

#### How to create a (good) new bug report:

Bugs are tracked as [GitHub issues](https://guides.github.com/features/issues/).

When you are creating a bug report, please do the following:

1. **Explain the problem**

   - *Use a clear and descriptive title* for the issue to identify the problem.
   - *Describe the exact steps which reproduce the problem* in as many details as possible.
     - Which function / method exactly you used?
     - What arguments or parameters were used?
     - *Provide a specific example*. (Includes links to pastebin, gists and so on.)
       If you're providing snippets in the issue, use
       [Markdown code blocks](https://help.github.com/articles/markdown-basics/#multiple-lines).

   - *Describe the behaviour you observed after following the steps*
     - Point out what exactly is the problem with that behaviour.
     - *Explain which behaviour you expected to see instead and why.*
     - *OPTIONALLY: Include screenshots and animated GIFs* which show you
       following the described steps and clearly demonstrate the problem.
       You can use [this tool](https://www.cockos.com/licecap/) to record GIFs on
       macOS and Windows, or [this tool](https://github.com/colinkeenan/silentcast)
       or [this tool](https://github.com/GNOME/byzanz) on Linux.

2. **Provide additional context for the problem (some of these may not always apply)**

   - *Did the problem start happening recently* (e.g. after updating to a new version)?
     - If the problem started recently, *can you reproduce the problem in older versions?*
     - Do you know the most recent package version in which the problem doesn't happen?

   - *Can you reliably reproduce the issue?* If not...
     - Provide details about how often the problem happens.
     - Provide details about under which conditions it normally happens.

   - Is the problem is related to *working with files*? If so....
     - Does the problem happen for all files and projects or only some?
     - Does the problem happen only when working with local or remote files?
     - Does the problem happen for files of a specific type, size, or encoding?
     - Is there anything else special about the files you are using?

3. **Include details about your configuration and environment**

- *Which version of the package are you using?*

- *What's the name and version of the OS you're using?*

- *Which julia packages do you have installed?*

- Are you using local configuration files to customize julia behaviour? If so...
  - Please provide the contents of those files, preferably in a
  [code block](https://help.github.com/articles/markdown-basics/#multiple-lines)
  or with a link to a [gist](https://gist.github.com/).


### Suggest an Enhancement

This section explains how to submit an enhancement proposal.
This includes completely new features, as well as minor improvements to existing functionality.
Following these suggestions will help maintainers and the community understand
your suggestion :pencil: and find related suggestions :mag_right:.

#### Before Submitting An Enhancement Proposal

1. **Perform a cursory issue search** to see if the enhancement has already been suggested.
2. If it has not, open a new issue as per the guidance below.
3. If it has...
  1. Add a comment to the existing issue instead of opening a new one.
  2. If it was closed, take the time to understand why this was so (it's ok to
     ask! :) ), and consider whether anything has changed that makes the reason
     outdated. If you can think of a convincing reason to reconsider the
     enhancement, feel free to open a new issue as per the guidance below.

#### How to submit a (good) new enhancement proposal

Enhancement proposals are tracked as
[GitHub issues](https://guides.github.com/features/issues/).

1. **Explain the enhancement**
   - *Use a clear and descriptive title* for the issue to identify the suggestion.
   - *Provide a step-by-step description of the suggested enhancement* in as many details as possible.
   - *Provide specific examples to demonstrate the steps*.
     Include copy/pasteable snippets which you use in those examples, as
     [Markdown code blocks](https://help.github.com/articles/markdown-basics/#multiple-lines).

   - If you want to change current behaviour...
     - Describe the *current* behaviour.
     - *Explain which behaviour you expected* to see instead and *why*.
     - *Will the proposed change alter APIs or existing exposed methods/types?*
       If so, this may cause dependency issues and breakages, so the maintainer
       will need to consider this when versioning the next release.

   - *OPTIONALLY: Include screenshots and animated GIFs*.
     You can use [this tool](https://www.cockos.com/licecap/) to record GIFs on
     macOS and Windows, and [this tool](https://github.com/colinkeenan/silentcast)
     or [this tool](https://github.com/GNOME/byzanz) on Linux.

2. **Provide additional context for the enhancement**

   - *Explain why this enhancement would be useful* to users and
     isn't something that can or should be implemented as a separate package.

   - *Do you know of other projects where this enhancement exists?*

3. **Include details about your configuration and environment**

   - Specify which *version of the package* you're using.

   - Specify the *name and version of the OS* you're using.

### Making Pull Requests

All julia packages can be developed locally.
For information on how to do this, see this section of the julia
[documentation](https://docs.julialang.org/en/stable/manual/packages/#Package-Development-1).

Before you start working on code, it is often a good idea to open an enhancement
[suggestion](#suggest-an-enhancement)

Once you decide to start working on code, the first thing you should do is make
yourself an account on [Github](https://github.com).
The chances are you already have one if you've done coding before and wanted to
make any scripts or software from a science project public.

The first step to contributing is to find the repository for the package you want to modify
(Microbiome.jl [is here](https://github.com/EcoJulia/Microbiome.jl),
BiobakeryUtils.jl [is here](https://github.com/EcoJulia/BiobakeryUtils.jl)).
Hit the 'Fork' button on the repositories page to create a forked copy of the
package for your own Github account. This is your blank slate to work on, and
will ensure your work and experiments won't hinder other users of the released
and stable package.

From there you can clone your fork of the package and work on it on your
machine using git.
Here's an example of cloning, assuming you already forked `Microbiome.jl`:

```sh
git clone https://github.com/<YOUR_GITHUB_USERNAME_HERE>/Microbiome.jl.git
```

Git will download or "clone" your fork and put it in a folder called
BioSequences.jl it creates in your current directory.

It is beyond the scope of this document to describe good git and github use in
more specific detail, as the folks at Git and GitHub have already done that wonderfully
on their own sites. If you have additional questions,
please feel free to start a [discussion on Microbiome.jl](https://github.com/EcoJulia/Microbiome.jl/discussions/new).

#### How to make (good) code contributions and new Pull-Requests

1. **In your code changes**

   - **Branch properly!**
     - If you are making a bug-fix, then you need to checkout your bug-fix branch
       from the last release tag.
     - If you are making a feature addition or other enhancement, checkout your
       branch from master.
     - See [here](#a-suggested-branching-model) for more information (or ask a package maintainer :smile:).

   - Follow the [julia style guide](https://docs.julialang.org/en/stable/manual/style-guide/).

   - Follow the [additional style suggestions](#additional-julia-code-style-suggestions).

   - Follow the [julia performance tips](https://docs.julialang.org/en/stable/manual/performance-tips/).

   - Update and add docstrings for new code, consistent with the [documentation styleguide](https://docs.julialang.org/en/stable/manual/documentation/).

   - Update information in the documentation located in the `docs/src/`
     folder of the package/repository if necessary.

   - Ensure that unit tests have been added which cover your code changes.

   - All changes should be compatible with the latest stable version of
     Julia.

   - Please comment liberally for complex pieces of internal code to facilitate comprehension.

2. **In your pull request**

   - *Describe* the changes in the pull request

   - Provide a *clear, simple, descriptive title*.

   - Do not include issue numbers in the PR title.

   - If you have implemented *new features* or behaviour
     - *Provide a description of the addition* in as many details as possible.
     - *Provide justification of the addition*.
     - *Provide a runnable example of use of your addition*. This lets reviewers
       and others try out the feature before it is merged or makes it's way to release.

   - If you have *changed current behaviour*...
     - *Describe the behaviour prior to you changes*
     - *Describe the behaviour after your changes* and justify why you have made the changes.
     - *Does your change alter APIs or existing exposed methods/types?*
       If so, this may cause dependency issues and breakages, so the maintainer
       will need to consider this when versioning the next release.
     - If you are implementing changes that are intended to increase performance, you
       should provide the results of a simple performance benchmark exercise
       demonstrating the improvement. Especially if the changes make code less legible.

#### Reviews and merging

You can open a pull request early on and push changes to it until it is ready,
or you can do all your editing locally and make a pull request only when it is
finished - it is up to you.

When your pull request is ready on Github,
mention one of the maintainers of the repo in a comment e.g. `@kescobo`
and ask them to review it.
You can also use Github's review feature.
They will review the code and documentation in the pull request,
and will assess it.

Your pull request will be accepted and merged if:

1. The dedicated package maintainers approve the pull request for merging.
2. The automated build system confirms that all unit tests pass without any issues.

It may also be that the reviewers or package maintainers will want to you to make
changes to your pull request before they will merge it.
Take the time to understand why any such request has been made,
and freely discuss it with the reviewers.
Feedback you receive should be constructive and considerate
(also see [here](#etiquette-and-conduct)).

## Styleguides

### Git Commit messages

* Use the present tense ("Add feature" not "Added feature").
* Use the imperative mood ("Move cursor to..." not "Moves cursor to...").
* Limit the first line to 72 characters or less.
* Reference issues and pull requests liberally after the first line.
* Consider starting the commit message with an applicable emoji:
    * :art: `:art:` when improving the format/structure of the code
    * :racehorse: `:racehorse:` when improving performance
    * :memo: `:memo:` when writing docs
    * :penguin: `:penguin:` when fixing something on Linux
    * :apple: `:apple:` when fixing something on macOS
    * :checkered_flag: `:checkered_flag:` when fixing something on Windows
    * :bug: `:bug:` when fixing a bug
    * :fire: `:fire:` when removing code or files
    * :green_heart: `:green_heart:` when fixing the CI build
    * :white_check_mark: `:white_check_mark:` when adding tests
    * :arrow_up: `:arrow_up:` when upgrading dependencies
    * :arrow_down: `:arrow_down:` when downgrading dependencies
    * :exclamation: `:exclamation:` when removing warnings or depreciations

### Additional julia style suggestions

- Indent with 4 spaces.

- For functions that are not a single expression, it is preferred to use an explicit `return`.
  Be aware that functions in julia implicitly return the the result of the last
  expression in the function, so plain `return` should be used to indicate that
  the function returns `nothing`.

- Type names are camel case, with the first letter capitalized. E.g.
  `SomeVeryUsefulType`.

- Module names should be camel case.

- Separate logical blocks of code with one blank line. Although it is common
  and acceptable for short single-line functions to be defined together on
  consecutive lines with no blank lines between them.

- Function names, apart from constructors, are all lowercase.
  Include underscores between words only if the name would be hard
  to read without.
  E.g.  `start`, `stop`, `find_letter` `find_last_digit`.
  It is good to separate concepts in a name with a `_`.

- Generally try to keep lines below 100-columns, unless splitting a long line
  onto multiple lines makes it harder to read.

- Files that declare modules should only declare the module, and import any
  modules that it requires. Any subsequent significant code should be included
  from separate files. E.g.

```julia
module AwesomeFeatures

using IntervalsTrees, JSON

include("feature1.jl")
include("feature2.jl")

end
```

- Files that declare modules should have the same name name of the module.
  E.g the module `SomeModule` is declared under the file `SomeModule.jl`.

- When extending method definitions, define the methods with a module name prefix. E.g.

```julia
function Base.start(iter::YourType)
  ...
end

Base.done(iter::YourType, state) = ...
```

- Functions that get or set variables in a struct should not be
  prefixed with 'get' or 'set'.
  The getter should be named for the variable it gets, and the setter
  should have the same name as the getter, with the suffix `!`.
  For example, for the variable `names`:

```julia
name(node) # get node name
name!(node, "somename") # set node name
```

- When using conditional branching, if code is statement-like, an
  if-else block should be used. However if the code is expression-like
  then julia's ternary operator should be used.
  ```julia
  matches == sketchlen ? 1.0 : matches / (2 * sketchlen - matches)
  ```
  Some simple checks and expressions are also expressed using the `&&` or `||`
  operators instead of if-else syntax. For example:
  ```julia
  isvalid(foo) || throw(ArgumentError("$foo is not valid"))
  ```
