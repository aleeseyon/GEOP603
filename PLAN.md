# GEOP 603 term project plan

Plan, work queue, run log and decisions for the GEOP 603 term project. This
project is a separate objective from any other research; nothing here feeds
into, or is driven by, other work.

## Ground rules
1. No facts from recall: every equation, constant and claim about prior work
   is backed by a citation verified in the session, with its DOI (digital
   object identifier) in `report/refs.bib`. Unverifiable items are marked
   `UNVERIFIED:` and listed in DECISIONS_NEEDED.
2. All numbers come from executed code; notebooks are executed end-to-end
   (`jupyter nbconvert --execute`) before commit.
3. Each new function gets a test against an exact or analytic limit first.
4. Define every acronym at first use in every document.
5. Surgical changes only; never delete existing code (comment out and flag).
6. Numerical QC at the end of every item: stability (CFL,
   Courant–Friedrichs–Lewy), dispersion, units, signs.
7. Public repository: no personal information.

## Project

Song & Fomel (2011), Fourier finite-difference (FFD) wave propagation, Geophysics
76(5), T123–T129. Deliverables due December: executed notebooks, report,
20-min talk, defense Q&A.
Note (2026-11-10): the author's FFD files (BP 2007 TTI vz/vx/eta/theta .bin,
grid.json, `ffd_tti_devito.py`, notebooks) live in a separate interactive chat
on TTI dispersion. That chat concluded (sources: Song & Fomel 2011, Song,
Fomel & Ying 2013 GJI 193(2) doi:10.1093/gji/ggt017, Koene et al. 2018 GJI
213) that FFD dispersion is medium- and path-dependent, so no single post-hoc
Fourier filter removes it; in-loop fixes are lowrank FFD (LFFD), a
higher-order FD part, or several reference media. Re-verify those citations
before reuse. Do not duplicate that work; C3 may build on it once files arrive.
- [ ] C0. If the author's notebooks are in GEOP603: re-execute, log any
      failures, do not change physics. If absent: skip to C2/C5 and note it.
- [ ] C1. Verification suite: homogeneous-medium FFD vs analytic 2D Green's
      function; dispersion vs 2nd/4th-order FD (finite difference) at equal
      Δt; stability limit scan vs reference-velocity choice.
- [ ] C2. Report skeleton `report/report.md` (Abstract, Theory with numbered
      equations from the paper, Implementation, Verification, Replication vs
      published figures, Extension, Conclusions); verified references.
- [ ] C3. Extension: **topic to be chosen by the author** (blocked). An
      earlier idea linked the extension to other research; it was removed on
      2026-11-11 so this project stays independent. See DECISIONS_NEEDED.
- [ ] C4. Slide outline + speaker notes for a 20-min talk (≈ 14 slides).
- [ ] C5. Defense Q&A bank: 30 likely examiner questions with short,
      referenced answers (stability, dispersion, reference velocity, TTI
      (tilted transverse isotropy), cost, limitations, extension).

Final interactive step (with the author): build the report as a Doc and the
talk as a Slides deck from the repo content.

## Run log
(one line per run: date, item, result, checks run, open issues)
- 2026-11-11 | setup | Plan moved here from a shared sprint file so that this course project is tracked separately | no code changed | Open: D1

## DECISIONS_NEEDED
- D1 (2026-11-11): choose the C3 extension topic (course-only).
