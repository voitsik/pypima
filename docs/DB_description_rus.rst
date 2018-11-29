Описание таблиц
===============

Таблицы с результатами подгонки лепестков
-----------------------------------------

Результаты подгонки лепестков хранятся в двух таблицах: ``pima_runs`` и ``pima_obs``. Таблица ``pima_runs`` содержит общую информацию об экспериментах, в частности время начала и окончания, параметры входных данных и т.д. Каждая запись в таблице ``pima_runs`` соответствует одной загрузке данных в **PIMA**, т.е. одному запуску таска ``load``. Большая часть данных для этой таблицы берётся из stt-файла.
Таблица ``pima_obs`` содержит непосредственно результаты подгонки лепестков. Каждая запись в таблице ``pima_obs`` соответствует одному "наблюдению" или одной строчке fri-файла, который является результатом работы таска ``frib``. Таблицы связаны через поля ``pima_runs.id`` и ``pima_obs.run_id`` так, что одной записи в ``pima_runs`` соответствуют несколько записей в ``pima_obs``.

Таблица ``pima_runs``
~~~~~~~~~~~~~~~~~~~~~

================== ===========================  ========================================
  Column                Type                     Описание
================== ===========================  ========================================
id                 integer                      Уникальный номер записи
exper_name         character varying(20)        Название эксперимента
band               band_type                    Код частотного диапазона (p, l, c, k)
fits_idi           character varying(256)       Имя FITS-IDI файла
scan_part          integer                      Номер "прогона" [1]_
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

.. [1] ``scan_part``: 1 - обработка полного скана, 2 - половины скана, 3 - треть и т.д.; 100 - длина скана 60 секунд.
   >100 - всякие тесты, для финальных результатов не используются.




Таблица ``pima_obs``
~~~~~~~~~~~~~~~~~~~~

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
st1          character varying(8)         Станция 1 (IVS Station name)
st2          character varying(8)         Станция 2 (IVS Station name)
delay        real                         Delay (s)
rate         real                         Rate (s/s)
accel        real                         Acceleration (s/s^2)
snr          real                         SNR
ampl         real                         Сырая амплитуда
solint       real                         Fringe fitting solution interval (s)
u            real                         U-координата (lambda)
v            real                         V-координата (lambda)
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


Таблица с результатами калибровки
---------------------------------

Результаты амплитудной калибровки записываются в таблицу ``ra_uvfits``. Эта таблица заполняется на основе UV-FITS файлов, полученных с помощью таска ``splt``.

Таблица ``ra_uvfits``
~~~~~~~~~~~~~~~~~~~~~

============ ============================= ==========================================
   Column                Type               Описание
============ ============================= ==========================================
 id          integer                        Уникальный номер записи
 source      character varying(20)          Имя источника (B1950)
 exper_name  character varying(20)          Название эксперимента
 band        band_type                      Код частотного диапазона (p, l, c, k)
 polar       polar_type                     Поляризация (RR, RL, LR, LL)
 sta1        character(2)                   Станция 1 (2-letter code)
 sta2        character(2)                   Станция 2 (2-letter code)
 u           real                           U-координата [2]_
 v           real                           V-координата [2]_
 ampl        real                           Калиброванная амплитуда (Jy)
 weight      real                           Вес ($1 / \sigma^2$)
 inttime     real                           Время интегрирования (s)
 file_name   character varying(256)         Имя UV-файла
 if_id       integer                        Номер частотного канала (IF)
 ind         integer                        Номер записи в UV-файле
 time        timestamp without time zone    Дата и время
 run_id      integer                        pima_runs.id
 freq        real                           Частота (Hz)
============ ============================= ==========================================

Примечания.

.. [2] Чтобы получить значения uv-координат в длинах волн, нужно домножить ``u`` и ``v`` на частоту ``freq``.


Вспомогательные таблицы
-----------------------
