import csv
import os
from itertools import islice
import datetime
import re


def main():
    data_root = os.path.abspath('../resource/raw_data/')
    mapping_root = os.path.abspath('../resource/reproduce_mapping/')
    save_root = os.path.abspath('../resource/reproduce_cache')

    diagnosis_path = os.path.join(data_root, 'diagnoses_icd.csv')
    admission_path = os.path.join(data_root, 'admissions.csv')
    lab_test_path = os.path.join(data_root, 'labevents.csv')
    medicine_path = os.path.join(data_root, 'prescriptions.csv')
    vital_sign_path = os.path.join(data_root, 'chartevents.csv')
    patient_path = os.path.join(data_root, 'patients.csv')
    operation_path = os.path.join(data_root, 'procedures_icd.csv')
    lab_mapping_path = os.path.join(mapping_root, 'lab_test_list.csv')
    diagnosis_mapping_path = os.path.join(mapping_root, 'disease_list.csv')
    medicine_mapping_path = os.path.join(mapping_root, 'drug_list.csv')
    operation_mapping_path = os.path.join(mapping_root, 'operation_list.csv')

    visit_dict = get_admissions(admission_path, save_root, read_from_cache=False)
    print('visit dict loaded')
    operation_dict = get_procedure(visit_dict, save_root, operation_path, operation_mapping_path, read_from_cache=False)
    print('operation dict loaded')
    event_dict = get_event(visit_dict, save_root, lab_test_path, read_from_cache=False)
    print('event dict loaded')
    age_sex_dict = get_sex_age(visit_dict, save_root, patient_path, read_from_cache=False)
    print('age sex dict loaded')
    vital_sign_dict = get_vital_sign(visit_dict, save_root, vital_sign_path, read_from_cache=False)
    print('vital sign dict loaded')
    medicine_dict = get_medicine(visit_dict, save_root, medicine_path, medicine_mapping_path, read_from_cache=False)
    print('medicine dict loaded')
    lab_test_dict = get_lab_test(visit_dict, save_root, lab_test_path, lab_mapping_path, read_from_cache=False)
    print('lab test dict loaded')
    diagnosis_dict = get_diagnosis(visit_dict, save_root, diagnosis_path, diagnosis_mapping_path, read_from_cache=False)
    print('diagnosis dict loaded')
    egfr_dict = egfr_calculation(visit_dict, age_sex_dict, lab_test_dict)
    print('egfr dict loaded')

    save_path = os.path.join(os.path.abspath('../resource/'), 'mimic_unpreprocessed.csv')
    reconstruct(visit_dict, operation_dict, event_dict, age_sex_dict, vital_sign_dict, medicine_dict, lab_test_dict,
                diagnosis_dict, egfr_dict, save_path)


def egfr_calculation(visit_dict, sex_age_dict, lab_test_dict):
    # based on the creatinine version of CKD-EPI equation
    egfr_dict = dict()
    for patient_id in visit_dict:
        egfr_dict[patient_id] = dict()
        for visit_id in visit_dict[patient_id]:
            egfr_dict[patient_id][visit_id] = -1

    for patient_id in egfr_dict:
        for visit_id in egfr_dict[patient_id]:
            ethnic_group = visit_dict[patient_id][visit_id]['ethnicity'].lower()
            try:
                creatinine = float(lab_test_dict[patient_id][visit_id]['SCr'][0])
            except ValueError:
                continue
            sex = float(sex_age_dict[patient_id][visit_id]['sex'])
            age = float(sex_age_dict[patient_id][visit_id]['age'])
            factor = 1
            if ethnic_group.__contains__('black') or ethnic_group.__contains__('africa'):
                factor = 1.159
            if creatinine <= 0 or sex < 0 or age <= 0:
                continue
            if sex == 0:
                if creatinine <= 0.7:
                    egfr_dict[patient_id][visit_id] = 144 * ((creatinine / 0.7) ** -0.329) * (0.993 ** age) * factor
                else:
                    egfr_dict[patient_id][visit_id] = 144 * ((creatinine / 0.7) ** -1.209) * (0.993 ** age) * factor
            elif sex == 1:
                if creatinine <= 0.9:
                    egfr_dict[patient_id][visit_id] = 141 * ((creatinine / 0.9) ** -0.411) * (0.993 ** age) * factor
                else:
                    egfr_dict[patient_id][visit_id] = 141 * ((creatinine / 0.9) ** -1.209) * (0.993 ** age) * factor
            else:
                raise ValueError('')
    return egfr_dict


def reconstruct(visit_dict, operation_dict, event_dict, age_sex_dict, vital_sign_dict, medicine_dict, lab_test_dict,
                diagnosis_dict, egfr_dict, save_path, one_visit=True):
    data_dict = dict()
    feature_list = ['egfr']
    for patient_id in visit_dict:
        data_dict[patient_id] = dict()
        visit_list = list()
        for visit_id in visit_dict[patient_id]:
            visit_list.append(int(visit_id))
        visit_list = sorted(visit_list)
        if one_visit:
            visit_list = [str(visit_list[-1])]
        else:
            visit_list = [str(item) for item in visit_list]

        for visit_id in visit_list:
            data_dict[patient_id][visit_id] = dict()
            data_dict[patient_id][visit_id]['egfr'] = egfr_dict[patient_id][visit_id]
            for operation in operation_dict[patient_id][visit_id]:
                data_dict[patient_id][visit_id][operation] = operation_dict[patient_id][visit_id][operation]
            for event in event_dict[patient_id][visit_id]:
                data_dict[patient_id][visit_id][event] = event_dict[patient_id][visit_id][event]
            for feature in age_sex_dict[patient_id][visit_id]:
                data_dict[patient_id][visit_id][feature] = age_sex_dict[patient_id][visit_id][feature]
            for vital_sign in vital_sign_dict[patient_id][visit_id]:
                data_dict[patient_id][visit_id][vital_sign] = vital_sign_dict[patient_id][visit_id][vital_sign]
            for medicine in medicine_dict[patient_id][visit_id]:
                data_dict[patient_id][visit_id][medicine] = medicine_dict[patient_id][visit_id][medicine]
            for lab_test in lab_test_dict[patient_id][visit_id]:
                data_dict[patient_id][visit_id][lab_test] = lab_test_dict[patient_id][visit_id][lab_test][0]
            for diagnosis in diagnosis_dict[patient_id][visit_id]:
                data_dict[patient_id][visit_id][diagnosis] = diagnosis_dict[patient_id][visit_id][diagnosis]

    for patient_id in visit_dict:
        for visit_id in visit_dict[patient_id]:
            for operation in operation_dict[patient_id][visit_id]:
                feature_list.append(operation)
            for event in event_dict[patient_id][visit_id]:
                feature_list.append(event)
            for feature in age_sex_dict[patient_id][visit_id]:
                feature_list.append(feature)
            for vital_sign in vital_sign_dict[patient_id][visit_id]:
                feature_list.append(vital_sign)
            for medicine in medicine_dict[patient_id][visit_id]:
                feature_list.append(medicine)
            for lab_test in lab_test_dict[patient_id][visit_id]:
                feature_list.append(lab_test)
            for diagnosis in diagnosis_dict[patient_id][visit_id]:
                feature_list.append(diagnosis)
            break
        break
    data_to_write = []
    head = ['patient_id', 'visit_id']
    for item in feature_list:
        head.append(item)
    data_to_write.append(head)
    for patient_id in data_dict:
        for visit_id in data_dict[patient_id]:
            line = [patient_id, visit_id]
            for index in range(2, len(head)):
                line.append(data_dict[patient_id][visit_id][head[index]])
            data_to_write.append(line)
    with open(save_path, 'w', encoding='utf-8-sig', newline='') as file:
        csv.writer(file).writerows(data_to_write)


def get_procedure(visit_dict, save_root, procedure_path, mapping_file, read_from_cache=True, file_name='procedure.csv'):
    if read_from_cache:
        procedure_dict = dict()
        with open(os.path.join(save_root, file_name), 'r', encoding='utf-8-sig', newline='') as file:
            csv_reader = csv.reader(file)
            for line in islice(csv_reader, 1, None):
                patient_id, visit_id, feature, value = line
                if not procedure_dict.__contains__(patient_id):
                    procedure_dict[patient_id] = dict()
                if not procedure_dict[patient_id].__contains__(visit_id):
                    procedure_dict[patient_id][visit_id] = dict()
                procedure_dict[patient_id][visit_id][feature] = int(value)
        return procedure_dict
    procedure_dict = dict()
    mapping_dict = dict()
    with open(mapping_file, 'r', encoding='utf-8-sig', newline='') as file:
        csv_reader = csv.reader(file)
        for line in islice(csv_reader, 1, None):
            name, code, _ = line
            mapping_dict[code] = name
    for patient_id in visit_dict:
        procedure_dict[patient_id] = dict()
        for visit_id in visit_dict[patient_id]:
            procedure_dict[patient_id][visit_id] = dict()
            for code in mapping_dict:
                procedure_dict[patient_id][visit_id][mapping_dict[code]] = 0

    with open(procedure_path, 'r', newline='') as file:
        csv_reader = csv.reader(file)
        for line in islice(csv_reader, 1, None):
            patient_id, visit_id, icd_9 = line[1], line[2], line[4]
            if not (procedure_dict.__contains__(patient_id) and procedure_dict[patient_id].__contains__(visit_id)
                    and mapping_dict.__contains__(icd_9)):
                continue
            procedure_dict[patient_id][visit_id][mapping_dict[icd_9]] = 1

    data_to_write = [['patient_id', 'visit_id', 'operation', 'positive']]
    for patient_id in procedure_dict:
        for visit_id in procedure_dict[patient_id]:
            for feature in procedure_dict[patient_id][visit_id]:
                value = procedure_dict[patient_id][visit_id][feature]
                data_to_write.append([patient_id, visit_id, feature, value])
    with open(os.path.join(save_root, file_name), 'w', encoding='utf-8-sig', newline='') as file:
        csv.writer(file).writerows(data_to_write)
    return procedure_dict


def get_event(visit_dict, save_root, lab_test_path, read_from_cache=True, file_name='event.csv'):
    if read_from_cache:
        event_dict = dict()
        with open(os.path.join(save_root, file_name), 'r', encoding='utf-8-sig', newline='') as file:
            csv_reader = csv.reader(file)
            for line in islice(csv_reader, 1, None):
                patient_id, visit_id, death, death_time, aki, aki_2_day, aki_7_day, aki_time = line
                if not event_dict.__contains__(patient_id):
                    event_dict[patient_id] = dict()
                event_dict[patient_id][visit_id] = {
                    'aki': int(aki), 'aki_time': float(aki_time), 'aki_2_day': int(aki_2_day),
                    'aki_7_day': int(aki_7_day), 'death': int(death), 'death_time': float(death_time)
                }
        return event_dict
    event_dict = dict()
    for patient_id in visit_dict:
        event_dict[patient_id] = dict()
        for visit_id in visit_dict[patient_id]:
            event_dict[patient_id][visit_id] = {'aki': 0,
                                                'aki_time': -1,
                                                'aki_2_day': 0,
                                                'aki_7_day': 0,
                                                'death': 0,
                                                'death_time': -1}

    # death evaluation
    for patient_id in visit_dict:
        for visit_id in visit_dict[patient_id]:
            admit_time = visit_dict[patient_id][visit_id]['admit_time']
            death_time = visit_dict[patient_id][visit_id]['death_time']
            if death_time < admit_time:
                continue
            else:
                time_diff = death_time-admit_time
                time_diff = (time_diff.days*24*3600 + time_diff.seconds)/24/3600
                event_dict[patient_id][visit_id]['death'] = 1
                event_dict[patient_id][visit_id]['death_time'] = time_diff

    # creatinine evaluation
    creatinine_dict = dict()
    with open(lab_test_path, 'r', newline='') as file:
        csv_reader = csv.reader(file)
        for line in islice(csv_reader, 1, None):
            patient_id, visit_id, lab_code, test_time, result = line[1: 6]
            if not (event_dict.__contains__(patient_id) and event_dict[patient_id].__contains__(visit_id)):
                continue
            # code of creatinine
            if lab_code != '50912' or len(test_time) < 10 or len(result) < 1:
                continue
            if not creatinine_dict.__contains__(patient_id):
                creatinine_dict[patient_id] = dict()
            if not creatinine_dict[patient_id].__contains__(visit_id):
                creatinine_dict[patient_id][visit_id] = list()
            test_time = datetime.datetime.strptime(test_time, '%Y-%m-%d %H:%M:%S')

            result_list = re.findall('[-+]?[\d]+(?:,\d\d\d)*[.]?\d*(?:[eE][-+]?\d+)?', result)
            if len(result_list) > 0:
                if result_list[0].__contains__(','):
                    result_list[0] = result_list[0].replace(',', '')
                result = float(result_list[0])
            else:
                continue
            creatinine_dict[patient_id][visit_id].append([float(result), test_time])
    for patient_id in creatinine_dict:
        for visit_id in creatinine_dict[patient_id]:
            creatinine_dict[patient_id][visit_id] = sorted(creatinine_dict[patient_id][visit_id], key=lambda x: x[1])

    for patient_id in creatinine_dict:
        for visit_id in creatinine_dict[patient_id]:
            creatinine_list = creatinine_dict[patient_id][visit_id]
            start_time = visit_dict[patient_id][visit_id]['admit_time']
            for index in range(len(creatinine_list)-1):
                value_1, test_time_1 = creatinine_list[index]
                if test_time_1 < start_time:
                    continue
                for index_ in range(index, len(creatinine_list)):
                    value_2, test_time_2 = creatinine_list[index_]
                    if value_2 == 0 or value_1 == 0:
                        continue
                    time_difference = (test_time_2 - test_time_1).days
                    if time_difference < 2:
                        if value_2 - value_1 > 0.3:
                            event_dict[patient_id][visit_id]['aki'] = 1
                            event_dict[patient_id][visit_id]['aki_2_day'] = 1
                            difference_from_start = ((test_time_2 - start_time).days * 3600 * 24 +
                                                     (test_time_2 - start_time).seconds) / (3600 * 24)
                            if event_dict[patient_id][visit_id]['aki_time'] == -1:
                                event_dict[patient_id][visit_id]['aki_time'] = difference_from_start
                            else:
                                if difference_from_start < event_dict[patient_id][visit_id]['aki_time']:
                                    event_dict[patient_id][visit_id]['aki_time'] = difference_from_start
                    if time_difference < 7:
                        if value_2/value_1 > 1.5:
                            event_dict[patient_id][visit_id]['aki'] = 1
                            event_dict[patient_id][visit_id]['aki_7_day'] = 1
                            difference_from_start = ((test_time_2 - start_time).days * 3600 * 24 +
                                                     (test_time_2 - start_time).seconds) / (3600 * 24)
                            if event_dict[patient_id][visit_id]['aki_time'] == -1:
                                event_dict[patient_id][visit_id]['aki_time'] = difference_from_start
                            else:
                                if difference_from_start < event_dict[patient_id][visit_id]['aki_time']:
                                    event_dict[patient_id][visit_id]['aki_time'] = difference_from_start

    data_to_write = [['patient_id', 'visit_id', 'death', 'death_time', 'aki', 'aki_2_day', 'aki_7_day', 'aki_time']]
    for patient_id in event_dict:
        for visit_id in event_dict[patient_id]:
            aki = event_dict[patient_id][visit_id]['aki']
            aki_time = event_dict[patient_id][visit_id]['aki_time']
            aki_2_day = event_dict[patient_id][visit_id]['aki_2_day']
            aki_7_day = event_dict[patient_id][visit_id]['aki_7_day']
            death = event_dict[patient_id][visit_id]['death']
            death_time = event_dict[patient_id][visit_id]['death_time']
            data_to_write.append([patient_id, visit_id, death, death_time, aki, aki_2_day, aki_7_day, aki_time])
    with open(os.path.join(save_root, file_name), 'w', encoding='utf-8-sig', newline='') as file:
        csv.writer(file).writerows(data_to_write)
    return event_dict


def get_sex_age(visit_dict, save_root, patient_path, read_from_cache=True, file_name='visit_info.csv'):
    if read_from_cache:
        sex_age_dict = dict()
        with open(os.path.join(save_root, file_name), 'r', newline='', encoding='utf-8-sig') as file:
            csv_reader = csv.reader(file)
            for line in islice(csv_reader, 1, None):
                patient_id, visit_id, sex, age = line
                if not sex_age_dict.__contains__(patient_id):
                    sex_age_dict[patient_id] = dict()
                sex_age_dict[patient_id][visit_id] = {'age': float(age), 'sex': int(sex)}
        return sex_age_dict
    sex_age_dict = dict()
    for patient_id in visit_dict:
        sex_age_dict[patient_id] = dict()
        for visit_id in visit_dict[patient_id]:
            sex_age_dict[patient_id][visit_id] = {'age': -1, 'sex': -1}

    with open(patient_path, 'r', newline='') as file:
        csv_reader = csv.reader(file)
        for line in islice(csv_reader, 1, None):
            patient_id, sex, birthday = line[1: 4]
            if not sex_age_dict.__contains__(patient_id):
                continue
            if len(sex) < 1 or len(birthday) < 10:
                continue
            birthday = datetime.datetime.strptime(birthday, '%Y-%m-%d %H:%M:%S')

            if sex == 'F':
                sex = 0
            elif sex == 'M':
                sex = 1
            else:
                raise ValueError('')
            for visit_id in visit_dict[patient_id]:
                admission_time = visit_dict[patient_id][visit_id]['admit_time']
                age = (admission_time-birthday).days / 365
                sex_age_dict[patient_id][visit_id] = {'age': age, 'sex': sex}

    data_to_write = [['patient_id', 'visit_id', 'age', 'sex']]
    with open(os.path.join(save_root, file_name), 'w', encoding='utf-8-sig', newline='') as file:
        for patient_id in sex_age_dict:
            for visit_id in sex_age_dict[patient_id]:
                data_to_write.append([patient_id, visit_id, sex_age_dict[patient_id][visit_id]['sex'],
                                      sex_age_dict[patient_id][visit_id]['age']])
        csv.writer(file).writerows(data_to_write)
    return sex_age_dict


def get_vital_sign(visit_dict, save_root, vital_sign_path, read_from_cache=True, file_name='vital_sign.csv'):
    if read_from_cache:
        vital_sign_dict = dict()
        with open(os.path.join(save_root, file_name), 'r', encoding='utf-8-sig', newline='') as file:
            csv_reader = csv.reader(file)
            for line in islice(csv_reader, 1, None):
                patient_id, visit_id, feature, value = line
                if not vital_sign_dict.__contains__(patient_id):
                    vital_sign_dict[patient_id] = dict()
                if not vital_sign_dict[patient_id].__contains__(visit_id):
                    vital_sign_dict[patient_id][visit_id] = dict()
                vital_sign_dict[patient_id][visit_id][feature] = float(value)
        return vital_sign_dict

    vital_sign_dict = dict()
    for patient_id in visit_dict:
        vital_sign_dict[patient_id] = dict()
        for visit_id in visit_dict[patient_id]:
            vital_sign_dict[patient_id][visit_id] = {
                'DBP': [-1, datetime.datetime(2500, 1, 1, 0, 0, 0, 0)],
                'SBP': [-1, datetime.datetime(2500, 1, 1, 0, 0, 0, 0)],
                'height': [-1, datetime.datetime(2500, 1, 1, 0, 0, 0, 0)],
                'weight': [-1, datetime.datetime(2500, 1, 1, 0, 0, 0, 0)],
            }

    with open(vital_sign_path, 'r', newline='', buffering=-1) as file:
        csv_reader = csv.reader(file)
        for line in islice(csv_reader, 1, None):
            patient_id, visit_id, item_id, chart_time, value, unit = \
                line[1], line[2], line[4], line[5], line[9], line[10]
            if not (vital_sign_dict.__contains__(patient_id) and vital_sign_dict[patient_id].__contains__(visit_id)):
                continue
            if len(chart_time) < 10 or len(value) < 1:
                continue
            value = float(value)
            unit = unit.lower()
            chart_time = datetime.datetime.strptime(chart_time, '%Y-%m-%d %H:%M:%S')
            # 1 lbs = 0.453592 kg
            # 1 inches = 2.54 cm
            # 1 feet = 30.48 cm
            # 1 oz = 0.0283495 kg
            # SBP
            if item_id == '51' or item_id == '455' or item_id == '220179' or item_id == '220050':
                if unit != 'mmhg':
                    continue
                if vital_sign_dict[patient_id][visit_id]['SBP'][1] > chart_time:
                    vital_sign_dict[patient_id][visit_id]['SBP'] = value, chart_time
            # dbp
            if item_id == '8368' or item_id == '8441' or item_id == '220180' or item_id == '220051':
                if unit != 'mmhg':
                    continue
                if vital_sign_dict[patient_id][visit_id]['DBP'][1] > chart_time:
                    vital_sign_dict[patient_id][visit_id]['DBP'] = value, chart_time
            # height
            if item_id == '216' or item_id == '1394' or item_id == '226707' or item_id == '226730' or item_id == '920':
                if unit == 'cm':
                    value = value
                elif unit == 'inch' or unit == 'inches':
                    value = value * 2.54
                elif unit == 'feet' or unit == 'feets':
                    value = value * 30.48
                else:
                    continue
                if vital_sign_dict[patient_id][visit_id]['height'][1] > chart_time and 250 > value > 50:
                    vital_sign_dict[patient_id][visit_id]['height'] = value, chart_time
            # weight
            if item_id == '3580' or item_id == '3581' or item_id == '3582' or item_id == '224639' or item_id == '763' \
                    or item_id == '226512' or item_id == '226531' or item_id == '762':
                if unit == 'kg':
                    value = value
                elif unit == 'lbs' or item_id == '226531':
                    value = value * 0.453592
                elif unit == 'oz':
                    value = value * 0.0283495
                else:
                    continue
                if vital_sign_dict[patient_id][visit_id]['weight'][1] > chart_time and 300 > value > 20:
                    vital_sign_dict[patient_id][visit_id]['weight'] = value, chart_time

    for patient_id in vital_sign_dict:
        for visit_id in vital_sign_dict[patient_id]:
            weight = vital_sign_dict[patient_id][visit_id]['weight'][0]
            height = vital_sign_dict[patient_id][visit_id]['height'][0]
            if weight == -1 or height == -1:
                bmi = -1
            else:
                bmi = weight * 10000 / height / height
            vital_sign_dict[patient_id][visit_id]['BMI'] = bmi, -1
            for feature in vital_sign_dict[patient_id][visit_id]:
                vital_sign_dict[patient_id][visit_id][feature] = vital_sign_dict[patient_id][visit_id][feature][0]

    data_to_write = [['patient_id', 'visit_id', 'feature', 'value']]
    with open(os.path.join(save_root, file_name), 'w', encoding='utf-8-sig', newline='') as file:
        for patient_id in vital_sign_dict:
            for visit_id in vital_sign_dict[patient_id]:
                for feature in vital_sign_dict[patient_id][visit_id]:
                    value = vital_sign_dict[patient_id][visit_id][feature]
                    data_to_write.append([patient_id, visit_id, feature, value])
        csv.writer(file).writerows(data_to_write)
    return vital_sign_dict


def get_medicine(visit_dict, save_root, medicine_path, mapping_file, read_from_cache=True, file_name='medicine.csv',
                 off_set=48):
    if read_from_cache:
        medicine_dict = dict()
        with open(os.path.join(save_root, file_name), 'r', encoding='utf-8-sig', newline='') as file:
            csv_reader = csv.reader(file)
            for line in islice(csv_reader, 1, None):
                patient_id, visit_id, drug, usage = line
                if not medicine_dict.__contains__(patient_id):
                    medicine_dict[patient_id] = dict()
                if not medicine_dict[patient_id].__contains__(visit_id):
                    medicine_dict[patient_id][visit_id] = dict()
                medicine_dict[patient_id][visit_id][drug] = int(usage)
        return medicine_dict

    medicine_dict = dict()
    name_cate_dict = dict()
    with open(mapping_file, 'r', encoding='utf-8-sig', newline='') as file:
        csv_reader = csv.reader(file)
        for line in csv_reader:
            category, _, english_name = line
            if name_cate_dict.__contains__(english_name.lower()):
                name_cate_dict[english_name.lower()].append(category)
            else:
                name_cate_dict[english_name.lower()] = [category]
    for patient_id in visit_dict:
        medicine_dict[patient_id] = dict()
        for visit_id in visit_dict[patient_id]:
            medicine_dict[patient_id][visit_id] = dict()
            for key in name_cate_dict:
                for item in name_cate_dict[key]:
                    medicine_dict[patient_id][visit_id][item] = 0

    with open(medicine_path, 'r', encoding='utf-8-sig', newline='') as file:
        csv_reader = csv.reader(file)
        for line in islice(csv_reader, 1, None):
            patient_id, visit_id, start_date = line[1], line[2], line[4]
            if len(start_date) < 10:
                continue
            if not (visit_dict.__contains__(patient_id) and visit_dict[patient_id].__contains__(visit_id)):
                continue
            start_date = datetime.datetime.strptime(start_date, '%Y-%m-%d %H:%M:%S')
            time_difference = (start_date - visit_dict[patient_id][visit_id]['admit_time'])
            if (time_difference.seconds+time_difference.days*24*3600)/3600 > off_set:
                continue

            drug_name = (line[7]+"_"+line[8]+'_'+line[9]).lower()
            for key in name_cate_dict:
                if drug_name.__contains__(key):
                    for item in name_cate_dict[key]:
                        medicine_dict[patient_id][visit_id][item] = 1

    data_to_write = [['patient_id', 'visit_id', 'medicine', 'usage']]
    for patient_id in medicine_dict:
        for visit_id in medicine_dict[patient_id]:
            for drug in medicine_dict[patient_id][visit_id]:
                data_to_write.append([patient_id, visit_id, drug, medicine_dict[patient_id][visit_id][drug]])
    with open(os.path.join(save_root, file_name), 'w', encoding='utf-8-sig', newline='') as file:
        csv.writer(file).writerows(data_to_write)
    return medicine_dict


def get_lab_test(visit_dict, save_root, lab_test_path, mapping_file, read_from_cache=True, file_name='lab_test.csv'):
    if read_from_cache:
        lab_test_dict = dict()
        with open(os.path.join(save_root, file_name), 'r', encoding='utf-8-sig', newline='') as file:
            csv_reader = csv.reader(file)
            for line in islice(csv_reader, 1, None):
                patient_id, visit_id, feature, value, record_time = line
                record_time = datetime.datetime.strptime(record_time, '%Y-%m-%d %H:%M:%S')
                result_list = re.findall('[-+]?[\d]+(?:,\d\d\d)*[.]?\d*(?:[eE][-+]?\d+)?', value)
                if len(result_list) > 0:
                    value = float(result_list[0])
                if not lab_test_dict.__contains__(patient_id):
                    lab_test_dict[patient_id] = dict()
                if not lab_test_dict[patient_id].__contains__(visit_id):
                    lab_test_dict[patient_id][visit_id] = dict()
                lab_test_dict[patient_id][visit_id][feature] = [value, record_time]
        return lab_test_dict

    # 读取mapping_file
    mapping_dict = dict()
    with open(mapping_file, 'r', encoding='utf-8-sig', newline='') as file:
        csv_reader = csv.reader(file)
        for line in islice(csv_reader, 1, None):
            lab_label, lab_chinese_name, lab_code = line[0: 3]
            mapping_dict[lab_code] = lab_label, lab_chinese_name
    # 构建模板
    lab_test_dict = dict()
    for patient_id in visit_dict:
        lab_test_dict[patient_id] = dict()
        for visit_id in visit_dict[patient_id]:
            lab_test_dict[patient_id][visit_id] = dict()
            for code in mapping_dict:
                lab_test_dict[patient_id][visit_id][mapping_dict[code][1]] = \
                    [-1, datetime.datetime(2500, 1, 1, 0, 0, 0, 0)]

    with open(lab_test_path, 'r', newline='') as file:
        csv_reader = csv.reader(file)
        for line in islice(csv_reader, 1, None):
            patient_id, visit_id, lab_code, test_time, result = line[1: 6]
            if not (lab_test_dict.__contains__(patient_id) and lab_test_dict[patient_id].__contains__(visit_id)):
                continue
            if (not mapping_dict.__contains__(lab_code)) or len(test_time) < 10:
                continue
            test_time = datetime.datetime.strptime(test_time, '%Y-%m-%d %H:%M:%S')

            if test_time < lab_test_dict[patient_id][visit_id][mapping_dict[lab_code][1]][1]:
                result_list = re.findall('[-+]?[\d]+(?:,\d\d\d)*[.]?\d*(?:[eE][-+]?\d+)?', result)
                if len(result_list) > 0:
                    if result_list[0].__contains__(','):
                        result_list[0] = result_list[0].replace(',', '')
                    result = float(result_list[0])
                lab_test_dict[patient_id][visit_id][mapping_dict[lab_code][1]] = [result, test_time]

    data_to_write = [['patient_id', 'visit_id', 'feature', 'value', 'record_time']]
    for patient_id in lab_test_dict:
        for visit_id in lab_test_dict[patient_id]:
            for feature in lab_test_dict[patient_id][visit_id]:
                value, record_time = lab_test_dict[patient_id][visit_id][feature]
                data_to_write.append([patient_id, visit_id, feature, value, record_time])
    with open(os.path.join(save_root, file_name), 'w', encoding='utf-8-sig', newline='') as file:
        csv.writer(file).writerows(data_to_write)
    return lab_test_dict


def get_admissions(admission_path, save_root, read_from_cache=True, file_name='admission.csv'):
    if read_from_cache:
        visit_dict = dict()
        with open(os.path.join(save_root, file_name), 'r', newline='', encoding='utf-8-sig') as file:
            csv_reader = csv.reader(file)
            for line in islice(csv_reader, 1, None):
                patient_id, visit_id, admit_time, discharge_time, death_time, ethnicity = line
                admit_time = datetime.datetime.strptime(admit_time, '%Y-%m-%d %H:%M:%S')
                discharge_time = datetime.datetime.strptime(discharge_time, '%Y-%m-%d %H:%M:%S')
                death_time = datetime.datetime.strptime(death_time, '%Y-%m-%d %H:%M:%S')
                if not visit_dict.__contains__(patient_id):
                    visit_dict[patient_id] = dict()
                visit_dict[patient_id][visit_id] = {'admit_time': admit_time, 'discharge_time': discharge_time,
                                                    'death_time': death_time, 'ethnicity': ethnicity}
        return visit_dict
    visit_dict = dict()
    with open(admission_path, 'r', newline='') as file:
        csv_reader = csv.reader(file)
        for line in islice(csv_reader, 1, None):
            patient_id, visit_id, admit_time, discharge_time, death_time = line[1: 6]
            ethnicity = line[13]

            admit_time = datetime.datetime.strptime(admit_time, '%Y-%m-%d %H:%M:%S')
            discharge_time = datetime.datetime.strptime(discharge_time, '%Y-%m-%d %H:%M:%S')
            if len(death_time) > 0:
                death_time = datetime.datetime.strptime(death_time, '%Y-%m-%d %H:%M:%S')
            else:
                death_time = datetime.datetime.strptime('1900-01-01 00:00:00', '%Y-%m-%d %H:%M:%S')
            if not visit_dict.__contains__(patient_id):
                visit_dict[patient_id] = dict()
            visit_dict[patient_id][visit_id] = {'admit_time': admit_time, 'discharge_time': discharge_time,
                                                "death_time": death_time, "ethnicity": ethnicity}

    with open(os.path.join(save_root, file_name), 'w', newline='', encoding='utf-8-sig') as file:
        data_to_write = [['patient_id', 'visit_id', 'admit_time', 'discharge_time', 'death_time', 'ethnicity']]
        for patient_id in visit_dict:
            for visit_id in visit_dict[patient_id]:
                admit_time = visit_dict[patient_id][visit_id]['admit_time']
                discharge_time = visit_dict[patient_id][visit_id]['discharge_time']
                death_time = visit_dict[patient_id][visit_id]['death_time']
                ethnicity = visit_dict[patient_id][visit_id]['ethnicity']
                data_to_write.append([patient_id, visit_id, admit_time, discharge_time, death_time, ethnicity])
        csv.writer(file).writerows(data_to_write)

    return visit_dict


def get_diagnosis(visit_dict, save_root, diagnosis_path, mapping_file, read_from_cache=True, file_name='diagnosis.csv'):
    if read_from_cache:
        diagnosis_dict = dict()
        with open(os.path.join(save_root, file_name), 'r', encoding='utf-8-sig', newline='') as file:
            csv_reader = csv.reader(file)
            for line in islice(csv_reader, 1, None):
                patient_id, visit_id, disease, positive = line
                if not diagnosis_dict.__contains__(patient_id):
                    diagnosis_dict[patient_id] = dict()
                if not diagnosis_dict[patient_id].__contains__(visit_id):
                    diagnosis_dict[patient_id][visit_id] = dict()
                diagnosis_dict[patient_id][visit_id][disease] = int(positive)
        return diagnosis_dict

    diagnosis_dict = dict()
    diagnosis_map_list = list()
    diagnosis_set = set()
    with open(mapping_file, 'r', encoding='utf-8-sig', newline='') as file:
        csv_reader = csv.reader(file)
        for line in islice(csv_reader, 1, None):
            diagnosis_map_list.append([line[1], line[4]])
            diagnosis_set.add(line[1])

    # 构建基本映射
    for patient_id in visit_dict:
        diagnosis_dict[patient_id] = dict()
        for visit_id in visit_dict[patient_id]:
            diagnosis_dict[patient_id][visit_id] = dict()
            for item in diagnosis_set:
                diagnosis_dict[patient_id][visit_id][item] = 0

    with open(diagnosis_path, 'r', newline='') as file:
        csv_reader = csv.reader(file)
        for line in islice(csv_reader, 1, None):
            _, patient_id, visit_id, _, icd_code = line
            if not (diagnosis_dict.__contains__(patient_id) and diagnosis_dict[patient_id].__contains__(visit_id)):
                continue
            for item in diagnosis_map_list:
                name, code = item
                if icd_code.__contains__(code):
                    diagnosis_dict[patient_id][visit_id][name] = 1

    data_to_write = [['patient_id', 'visit_id', 'disease', 'positive']]
    for patient_id in diagnosis_dict:
        for visit_id in diagnosis_dict[patient_id]:
            for disease in diagnosis_dict[patient_id][visit_id]:
                data_to_write.append([patient_id, visit_id, disease, diagnosis_dict[patient_id][visit_id][disease]])
    with open(os.path.join(save_root, file_name), 'w', encoding='utf-8-sig', newline='') as file:
        csv.writer(file).writerows(data_to_write)
    return diagnosis_dict


if __name__ == '__main__':
    main()
