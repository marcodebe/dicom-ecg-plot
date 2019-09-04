# -*- coding: utf-8 -*-
import os
import locale
import gettext

APP_NAME = "ecg"
LOCAL_PATH = 'locale/'

def get_lang():
    langs = []
    lc, encoding = locale.getdefaultlocale()

    if (lc):
        langs = [lc]

    language = os.environ.get('LANGUAGE', None)
    if (language):
        langs += language.split(":")

    langs += ["it_IT", "en_US"]

    return langs

gettext.bindtextdomain(APP_NAME, LOCAL_PATH)
gettext.textdomain(APP_NAME)

lang = gettext.translation(APP_NAME, LOCAL_PATH,
                           languages=get_lang(), fallback=True)

_ = lang.gettext

# Down here the translated words
ventr_freq = _('Ventr. Freq.')
pr_interval = _('PR Interval')
qrs_duration = _('QRS Duration')
qt_qtc = _('QT/QTc')
prt_axis = _('P-R-T Axis')
pat_id = _('Pat. ID')
pat_sex = _('sex')
pat_bdate = _('birthdate')
pat_age = _('year old')
duration = _('total time')
sampling_frequency = _('sample freq.')
acquisition_date = _('acquisition date')
