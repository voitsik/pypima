Описание таблиц
===============

Таблицы с результатами подгонки лепестков
-----------------------------------------

Table "public.pima_runs"

================== ===========================  ========================================
  Column                Type                     Описание
================== ===========================  ========================================
id                 integer                      Уникальный номер записи
exper_name         character varying(20)        Название эксперимента
band               band_type                    Код частотного диапазона (p, l, c, k)
fits_idi           character varying(256)       Имя FITS-IDI файла
scan_part          integer                      Номер "прогона" [1]
sp_chann_num       integer                      Количество спектральных каналов
time_epochs_num    integer                      Количество отсчетов по времени
scans_num          integer                      Количество сканов
obs_num            integer                      Количество "наблюдений"
uv_points_num      integer                      Общее количество записей в UV-данных
uv_points_used_num integer                      Количество используемых UV-записей
deselected_points  integer                      Общее количество исключенных записей
no_auto_points_num integer                      Количество записей, исключенных из-за
                                                отсутствия автоспектра
accum_length       real                         Время накопления в корреляторе
utc_minus_tai      interval                     UTC-TAI в секундах
nominal_start      timestamp without time zone  Дата и время начала эксперимента
nominal_end        timestamp without time zone  Дата и время окончания эксперимента
proc_date          timestamp without time zone  Дата и время обработки
last_error         character varying(256)       Информация об ошибке обработки
hostname           character varying(64)        Имя хоста на котором проводилась
                                                обработка
pima_version       character varying(8)         Версия PIMA
correlator_name    character varying(8)         Название коррелятора (обычно ASCFX)
================== ===========================  ========================================

Примечания.

.. [1] `scan_part`: 1 - обработка полного скана, 2 - половины скана, 3 - треть и т.д.; 100 - длина скана 60 секунд.
   >100 - всякие тесты, для финальных результатов не используются.




Table "public.pima_obs"

===========  =========================== ========================================
   Column               Type              Описание
===========  =========================== ========================================
id           integer                      Уникальный номер записи
obs          smallint                     Номер "наблюдения" в PIMA
scan_name    character varying(10)        Код скана (ddd-hhmm)
start_time   timestamp without time zone  Дата и время начала наблюдения
stop_time    timestamp without time zone  Дата и время конца наблюдения
source       character varying(20)        Имя источника (IVS)
polar        polar_type                   Поляризация
st1          character varying(8)         Станция 1
st2          character varying(8)         Станция 2
delay        real                         Delay (s)
rate         real                         Rate (s/s)
accel        real                         Acceleration (s/s^2)
snr          real                         SNR
ampl         real                         Сырая амплитуда
solint       real                         Fringe fitting solution interval (s)
u            real                         U
v            real                         V
base_ed      real                         Проекция базы в диаметрах Земли
ref_freq     real                         Частота наблюдения (Hz)
exper_name   character varying(20)        Название эксперимента
band         band_type                    Код частотного диапазона (p, l, c, k)
status       status_type                  Статус детектирования (y, u, n)
run_id       integer                      pima_runs.id
if_id        smallint                     Номер IF -- должно быть 0, поскольку используются все IF сразу
elevation    real[]                       Высота источника над горизонтом (массив из двух элементов)
bandpass     boolean                      Флаг bandpass калибровки
===========  =========================== ========================================
