<h1 align="center">
  <br>
  [Project's Logo]
  <br>
  QUICK-QM
  <br>
</h1>

<h4 align="center">This Python based project serves as a wrapper around the ORCA tool suite to quickly retrieve molecular properties and interaction energies of two molecular systems.</h4>

<!--<p align="center">
  [Project's badges]
</p>-->

<p align="center">
  <a href="#about">About</a> •
  <a href="#key-features">Key Features</a> •
  <a href="#getting-started">Getting Started</a> •
  <a href="#faq">FAQ</a> •
  <a href="#roadmap">Roadmap</a> •
  <a href="#support">Support</a> •
  <a href="#license">License</a>
</p>

![screenshot](screenshots/1.jpg)

## About

The ORCA tool suite offers various functions to either retrieve properties of a single molecular (A) or retrieving the interaction energy forming an adduct of n molecules (B), which both are commonly used approaches in various QM based studies.

While (A) requires a simple run and adequate calling of sub procedures, (B) requires first the calculation of the individual system energies and then retrieval of the interaction energy by forming an adduct of n molecules and running a calculation afterwards. The final interaction energy is given by the difference of the adduct system energy and the sum of all single system energies.

As (A) is a key component of any single molecular property retrieval and (B) could be used to model any molecule-molecule interaction, this approach can overall be widely used.

This package provides an easy-to-use wrapper, especially targeting two typical conducting simulation approaches: Single-Molecule-QM, to retrieve properties of a single molecule and Two-Molecule-QM, allowing to simulate the interaction of two molecules.

The tool-suite reduces the complexity and speeds up the process by defining the ORCA input files and allows a selector based alignment, quickly creating the alignments needed based on molecular functional groups. 

## Key Features

- Feature 1
- Feature 2
  - Feature 2.1
  - Feature 2.2
- Feature 3

## Getting Started

This is an example of how you may give instructions on setting up your project locally. To get a local copy up and running follow these simple example steps.

### Prerequisites

This is an example of how to list things you need to use the software and how to install them.

- npm

### Installing and Running

Below is an example of how you can instruct your audience on installing and setting up your app. This template doesn't rely on any external dependencies or services.

```bash
# Clone this repository
$ git clone https://github.com/[your-github-user]/[your-repository-name]

# Go into the repository
$ cd [your-repository-name]
```

## FAQ

### Is it any good?

[yes.](https://news.ycombinator.com/item?id=3067434)

## Definitions

- This are constraints you give to your project
- Here's another constraint: use this README!

## Roadmap

- [x] Add this README.
- [ ] Finish the unicorn project.

## Emailware

<!--[Project's name] is an [emailware](https://en.wiktionary.org/wiki/emailware). Meaning, if you liked using this app or it has helped you in any way, I'd like you send me an email at <your-email@email.com> about anything you'd want to say about this software. I'd really appreciate it!-->

## Contributing

If you have a suggestion that would make this better, please fork the repo and create a pull request. You can also simply open an issue with the tag "enhancement".
Don't forget to give the project a star! Thanks again!

1. Fork the project
2. Clone it open the repository in command line
3. Create your feature branch (`git checkout -b feature/amazing-feature`)
4. Commit your changes (`git commit -m 'Add some amazing feature'`)
5. Push to the branch (`git push origin feature/amazing-feature`)
6. Open a pull request from your feature branch from your repository into this repository main branch, and provide a description of your changes

## Support

You can also support us by:

<p align="left">
  <a href="https://buymeacoffee.com/bastian.eisenmann" target="_blank"><img src="https://www.buymeacoffee.com/assets/img/custom_images/purple_img.png" alt="Buy Me A Coffee" style="height: 41px !important;width: 174px !important;box-shadow: 0px 3px 2px 0px rgba(190, 190, 190, 0.5) !important;-webkit-box-shadow: 0px 3px 2px 0px rgba(190, 190, 190, 0.5) !important;" ></a> &nbsp &nbsp
  <!--<a href="https://www.patreon.com">
    <img src="https://c5.patreon.com/external/logo/become_a_patron_button@2x.png" width="160">
  </a>-->
</p>

## License (BSD 3-Clause License)

Copyright 2026 Eisenmann, Bastian

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS “AS IS” AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

<!--
## Acknowledgments

Inspiration, code snippets, etc.

- [Markdownify's README](https://github.com/amitmerchant1990/electron-markdownify#readme)

## You may also like...

List of apps or libs that do similar stuff as your project.

- [Best-README-Template](https://github.com/othneildrew/Best-README-Template)
- [Simple README.md template](https://gist.github.com/DomPizzie/7a5ff55ffa9081f2de27c315f5018afc)

---
-->

> [hidek.xyz](https://hidek.xyz) &nbsp;&middot;&nbsp;
> GitHub [@abehidek](https://github.com/abehidek) &nbsp;&middot;&nbsp;
> Twitter [@guilhermehabe](https://twitter.com/guilhermehabe)
