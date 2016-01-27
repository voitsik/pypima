2.3 (unreleased)
----------------


2.2 (2016-01-27)
----------------

New Features
^^^^^^^^^^^^

- Use ``logging`` module. Add `--log-file` command-line option to ``pima_auto.py``
  script.

- Add `--force-small` command-line option to ``pima_auto.py`` script to select
  FITS files correlated with small window only.

- Support new ASC FTP archive structure.

Other Changes and Additions
^^^^^^^^^^^^^^^^^^^^^^^^^^^

- Increase SNR detection limit to 5.9 for K-band.

- Comment out GAIN lines in ANTAB files for different frequency.
