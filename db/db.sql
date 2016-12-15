-- psql ra_results

CREATE TABLE vex_files (
    file_name varchar(32),
    exper_name varchar(32) primary key,
    exper_description varchar(128),
    PI_name varchar(64),
    PI_email varchar(64),
    exper_nominal_start timestamp,
    exper_nominal_stop timestamp
);
GRANT SELECT ON vex_files TO guest;

CREATE TABLE sources (
  IVS_name varchar(20) primary key,
  J2000_name varchar(20),
  B1950_name varchar(20),
  aips_name varchar(20),
  coordinates point not null,
  comments text,
  z real DEFAULT 0.0
);
GRANT SELECT ON sources TO guest;

CREATE TABLE pima_runs (
  id SERIAL PRIMARY KEY,
  exper_name varchar(20) references vex_files(exper_name),
  band char(1) NOT NULL,
  fits_idi varchar(256) NOT NULL,
  scan_part int NOT NULL,
  sp_chann_num int,
  time_epochs_num int,
  scans_num int,
  obs_num int,
  uv_points_num int,
  uv_points_used_num int,
  deselected_points_num int,
  no_auto_points_num int,
  accum_length real,
  utc_minus_tai interval,
  nominal_start timestamp,
  nominal_end timestamp,
  proc_date timestamp,
  last_error varchar(256),
  hostname varchar(64),
  pima_version varchar(8),
  UNIQUE (exper_name, band, fits_idi, scan_part)
);

CREATE INDEX pima_runs_exper_name_band_idx ON pima_runs (exper_name, band);

GRANT SELECT ON pima_runs TO guest;
GRANT SELECT, UPDATE, INSERT, DELETE ON pima_runs TO editor;
GRANT USAGE, SELECT ON SEQUENCE pima_runs_id_seq TO editor;

--
-- PIMA observations and fringe fitting results
--
CREATE TABLE pima_obs (
  id SERIAL primary key,
  obs smallint,
  scan_name varchar(10),
  start_time timestamp NOT NULL,
  stop_time timestamp NOT NULL,
  source varchar(20) references sources(IVS_name),
  polar char(2) NOT NULL,
  st1 varchar(8) NOT NULL,
  st2 varchar(8) NOT NULL,
  delay real,
  rate  real,
  accel real,
  snr real,
  ampl real,
  solint real,
  u real,
  v real,
  base_ed real,
  ref_freq real,
  exper_name varchar(20) NOT NULL,
  band char(1) NOT NULL,
  status char(1) DEFAULT 'u',
  run_id int references pima_runs(id) ON DELETE CASCADE,
  if_id smallint DEFAULT 0,
  elevation real[] DEFAULT ARRAY[0.0, 0.0],
  bandpass boolean
);

CREATE INDEX pima_obs_exper_name_band_idx ON pima_obs (exper_name, band);

GRANT SELECT ON pima_obs TO guest;
GRANT SELECT, UPDATE, INSERT, DELETE ON pima_obs TO editor;
GRANT USAGE, SELECT ON SEQUENCE pima_obs_id_seq TO editor;


CREATE TABLE clock_models (
    id SERIAL primary key,
    sta varchar(8) NOT NULL,
    time timestamp NOT NULL,
    clock_offset real,
    clock_rate real,
    group_delay real,
    delay_rate real,
    run_id int references pima_runs(id) ON DELETE CASCADE
);

GRANT SELECT ON clock_models TO guest;
GRANT SELECT, UPDATE, INSERT, DELETE ON clock_models TO editor;
GRANT USAGE, SELECT ON SEQUENCE clock_models_id_seq TO editor;

--
-- FITS-IDI files in ASC correlator archive
--
CREATE TABLE fits_files (
    path varchar(512) PRIMARY KEY,
    basename varchar(256) NOT NULL,
    exper_name varchar(32) REFERENCES vex_files(exper_name),
    obs_date date,
    corr_date date,
    correlator varchar(8),
    oper varchar(32),
    band char(1),
    ch_num int
    mdate date DEFAULT '2000-01-01',
    size bigint DEFAULT 0,
    comment varchar(64),
    ftp_user varchar(16)
);
CREATE INDEX fits_files_exper_name_band_idx ON fits_files (exper_name, band);

GRANT SELECT ON fits_files TO guest;
GRANT SELECT, UPDATE, INSERT, DELETE ON fits_files TO editor;

--
-- Calibrated UV-FITS files after PIMA processing
--
CREATE TABLE ra_uvfits (
    id SERIAL primary key,
    source varchar(20),
    exper_name varchar(20) references vex_files(exper_name),
    band char(1),
    polar char(2),
    sta1 char(2),
    sta2 char(2),
    u real,
    v real,
    freq real,
    ampl real,
    weight real,
    inttime real,
    file_name varchar(256),
    if_id int,
    ind int,
    time timestamp,
    run_id int references pima_runs(id) ON DELETE CASCADE
);

GRANT SELECT ON ra_uvfits TO guest;
GRANT SELECT, UPDATE, INSERT, DELETE ON ra_uvfits TO editor;
GRANT USAGE, SELECT ON SEQUENCE ra_uvfits_id_seq TO editor;

--
-- PIMA station names catalog
--
CREATE TABLE station_names (
    corr_name varchar(2) primary key,
    ivs_name varchar(8) NOT NULL,
    x real,
    y real,
    z real
);
GRANT SELECT ON station_names TO guest;
GRANT SELECT, UPDATE, INSERT, DELETE ON station_names TO editor;


--
-- Autocorrelation spectra
--
CREATE TABLE autospec_info (
    id serial primary key,
    exper_name varchar(20) references vex_files(exper_name),
    band char(1),
    polar char(2),
    sta varchar(8) NOT NULL,
    start_date timestamp,
    stop_date timestamp,
    run_id int references pima_runs(id) ON DELETE CASCADE
);
GRANT SELECT ON autospec_info TO guest;
GRANT SELECT, UPDATE, INSERT, DELETE ON autospec_info TO editor;
GRANT USAGE, SELECT ON SEQUENCE autospec_info_id_seq TO editor;

CREATE TABLE autospec (
    id bigserial primary key,
    if_num smallint,
    chann_num smallint,
    freq real,
    ampl real,
    info_id integer references autospec_info(id) ON DELETE CASCADE
);
GRANT SELECT ON autospec TO guest;
GRANT SELECT, UPDATE, INSERT, DELETE ON autospec TO editor;
GRANT USAGE, SELECT ON SEQUENCE autospec_id_seq TO editor;
