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
  UNIQUE (exper_name, band, fits_idi, scan_part)
);

CREATE INDEX pima_runs_exper_name_band_idx ON pima_runs (exper_name, band);

GRANT SELECT ON pima_runs TO guest;
GRANT SELECT, UPDATE, INSERT, DELETE ON pima_runs TO editor;
GRANT USAGE, SELECT ON SEQUENCE pima_runs_id_seq TO editor;

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
  run_id int references pima_runs(id) ON DELETE CASCADE
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

