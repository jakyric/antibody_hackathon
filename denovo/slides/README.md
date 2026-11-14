# Pipeline overview deck

Private Slides artifact (14 slides), created by the publish step; the URL is not stored here. Share it from the page's Share menu if others should open it.

Every number is read from `results/denovo/*` and the stage folders when the deck is built, so a refresh is one command plus one re-publish.

## Refresh

```bash
denovo/.venv/bin/python denovo/slides/build_slides.py --plan     # rewrites denovo/slides/artifact/project/slides/*.html, prints the publish plan
```

Then publish the slide files to the SAME artifact (Artifact tool, publish): `url` = the URL above, `root` = `denovo/slides/artifact` (absolute path),
`file_path` = `<root>/project/slides/cover.html`, `files` = the other thirteen `project/slides/<id>.html` entries listed in
`denovo/slides/publish_plan.json`. `project/deck.json` is part of the plan because the slide order changed (new `pose` slide); `project/deck.css` is only sent when the look changes.

Slide ids in order (final state, 14 slides): cover, task, pipeline, generate, validate, gate, scoring, pose, ranking, live, dist, final, limits, value.
`pose` compares the pose recovery of the scorers on the full HER2 domain and on the cropped HER2 539-607 re-run (reads `pose_check_vm.csv` and `crop/*.csv`, "pending" when the crop files are missing). `dist` is the boxplot figure (drawn from the CSVs as SVG, with the final 20 marked); `final` is the final list with its flags and the open question.

## Files

* `build_slides.py`: reads the result files, writes the slide HTML and `deck.json` (all content and charts live here).
* `artifact/project/deck.css`: the shared look (static). `artifact/project/slides/*.html` and `deck.json`: generated.
* `preview.py`: optional local layout check (headless Chromium, four slides per image); needs a working `chromium`.
* `publish_plan.json`: the exact publish parameters, rewritten on every build.

## What changes on refresh

Counts (designs, gate, scored, Tm, yield), the gate chart, the OpenFold3 and Tm strips, the candidate table, the Pareto scatter
(from `leaderboard_pareto.csv` and `submission_top20_pareto.csv`), JapanFold spend (from `status.md`), the composite weights
(parsed from `rank.py`). Anything without evidence yet shows "pending". No credentials or access details are read or written.
