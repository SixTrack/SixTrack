# Editing the SixTrack Manual

For version 5, the entire manual has been reformatted both for consistency, and to work with scripts that converts the LaTeX source to HTML.

This file describes some key points when adding to or editing the manual in the future.

## Source Formatting

* The LaTeX source is indented with 4 spaces. Please keep it consistent
* Use a line break after each sentence. This makes it much easier to track changes in git.
* Please add a non-breaking space `~` in references to avoid splitting between label and number. For example `Table~\ref{sometable}`.

## General Formatting Guidelines

* Avoid wrapping text in formats as `\emph{}`, `{\em }`, `\verb| |`, etc. Use instead `\textit{}`, `\textbf{}`, `\texttt{}`, etc. These format correctly, and among other things avoids the need for hacks like `\/` after italics text to retain proper word separation.
* *All* references to SixTrack keywords, to internal variables, and to file names, should be wrapped in `\texttt{}`.
* Hyphenation should use the proper hyphen character `-` (single -). Words like `order-of-magnitude` and `tune-shift` are hyphenated.
* Dashes should use the proper dash character `--` (double -). These are used to describe intervals, `5--10`, or to combine words that describes interactions, like `beam--beam`. It can also be used as a minus sign for negative numbers, but ideally those should be wrapped in `$$`.

## Keywords and Format Description

* The updated format of keyword descriptions uses the tabular environment rather than the previous method of paragraphs and subparagraphs. The first column should be wrapped in `\textbf{}`, and the description in regular font unless they describe keywords, in which case they should be rwapped in `\texttt{}`.

## Long Tables and Lists

* For tables that will use most of, or more than, a single page, use the `\begin{longtable}` environment instead. Please look at existing tables for how to format these. They differ from regular tables.
* For long Format Description section, there is also the option to use `\begin{longtabu}`, which automatically splits between pages which regular tabular does not.

## Structure

* Each new block should be described within a single `\section{}` in the appropriate chapter. Please add a `\label{}` as well.
* Long sections should be split into `\subsection{}` where appropriate.
* Do **not** use `\subsubsection{}`. These do not format well when converting to HTML. If you need headings in the text, one option is to use `\paragraph{}~\\`.
