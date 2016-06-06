---
layout: page
title: Home
permalink: /
---

This is the central repository for the MOSFIRE DRP developed by N. Konidaris and C. Steidel at Caltech, and currently hosted at the WMK Observatory.

If you need help with the pipeline or to report a problem, please visit our [issue tracking page](https://github.com/Keck-DataReductionPipelines/MosfireDRP/issues) hosted at GitHub. Please note that you need a free GitHub account to submit a ticket.

The currently release installation and reduction instructions are provided in the [DRP manual](manual).

The development and support team includes:

* Marc Kassis, Luca Rizzi, Jim Lyke, and Josh Walawender at W. M. Keck Observatory
* Nick Konidaris, Chuck Steidel at Caltech
* Tuan Do at UCLA

For direct communication with the support and development team, please email [mosfiredrp@gmail.com](mailto:mosfiredrp@gmail.com)

<p class="rss-subscribe">Subscribe to blog posts <a href="{{ "/feed.xml" | prepend: site.baseurl }}">via RSS</a></p>

<hr>

# MOSFIRE DRP Blog Posts

  <ul class="post-list">
    {% for post in site.posts %}
      <li>
        <span class="post-meta">{{ post.date | date: "%b %-d, %Y" }}</span>

        <h2>
          <a class="post-link" href="{{ post.url | prepend: site.baseurl }}">{{ post.title }}</a>
        </h2>
      </li>
    {% endfor %}
  </ul>

