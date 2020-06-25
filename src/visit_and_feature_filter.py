import os
import csv
import re
from itertools import islice


def main():
    patient_delete_missing_rate = 0.3
    egfr_threshold = 60
    item_list = ['egfr', 'age', 'sex', 'Angiography', 'PCI', 'diuretic', 'Anticoagulants', 'CCB', 'ACEI',
                 'PositiveInotropicDrugs', 'ARB', 'beta-blocker', 'antiplatelet', 'vasodilator', 'DBP', 'SBP',
                 'BMI', 'Triglycerides', 'TotalProtein', 'ALT', 'Sodium', 'GGT', 'Potassium',
                 'Glucose', 'TotalBilirubin', 'HDL-C', 'TnT', 'Hemoglobin', 'Calcium', 'NT-pro-BNP', 'Urea',
                 'AST', 'LDL-C', 'LVEF', 'VHD', 'myocardiopathy', 'CHD', 'stroke', 'AF']
    reserve_set = set(item_list)

    source_file_path = os.path.abspath('../resource/mimic_unpreprocessed.csv')
    target_file_path = os.path.join('../resource/', 'filtered.csv')
    index_item_dict = get_item_index_dict(source_file_path)
    feature_dict = read_un_preprocessed_data(source_file_path, index_item_dict)
    print('un preprocessed data size:' + str(calculate_visit_count(feature_dict)))
    feature_dict = unit_transform(feature_dict)
    print('unit transformed')
    feature_dict = discard_illegal_data_value(feature_dict)
    print('after delete illegal data value:' + str(calculate_visit_count(feature_dict)))
    feature_dict = delete_by_kidney_function(feature_dict, egfr_threshold)
    print('after delete by kidney data size:' + str(calculate_visit_count(feature_dict)))
    feature_dict = delete_by_admission_reason(feature_dict)
    print('after delete by admission reason, size:' + str(calculate_visit_count(feature_dict)))
    feature_dict = delete_juveniles(feature_dict)
    print('after delete juvenile data size:' + str(calculate_visit_count(feature_dict)))
    feature_dict = reserve_selected_feature(feature_dict, reserve_set)
    feature_dict = delete_visit_missing_too_much(feature_dict, patient_delete_missing_rate)
    print('after delete visit missing too much data, size:' + str(calculate_visit_count(feature_dict)))

    print('final size:' + str(calculate_visit_count(feature_dict)))

    head = ['patient_id', 'visit_id']
    for item_name in item_list:
        head.append(item_name)
    data_to_write = [head]
    for patient_id in feature_dict:
        for visit_id in feature_dict[patient_id]:
            line = list()
            info_dict = feature_dict[patient_id][visit_id]
            line.append(patient_id)
            line.append(visit_id)
            for item_name in item_list:
                if not info_dict.__contains__(item_name):
                    line.append(-1)
                else:
                    line.append(info_dict[item_name])
            data_to_write.append(line)
    with open(target_file_path, 'w', encoding='utf-8-sig', newline='') as file:
        csv.writer(file).writerows(data_to_write)


def reserve_selected_feature(feature_dict, reserve_set):
    new_feature_dict = dict()
    for patient_id in feature_dict:
        new_feature_dict[patient_id] = dict()
        for visit_id in feature_dict[patient_id]:
            new_feature_dict[patient_id][visit_id] = dict()
            for item in reserve_set:
                if feature_dict[patient_id][visit_id].__contains__(item):
                    new_feature_dict[patient_id][visit_id][item] = feature_dict[patient_id][visit_id][item]
    return new_feature_dict


def unit_transform(feature_dict):
    # Scr     1 mg/dL = 88.41 umol/L
    # Triglycerides  1 mmol/L = 88.6 mg/dl
    # Glucose 1 mmol/L = 18 mg/dL
    # TotalBilirubin 1 mg/dL = 17.1 umol/L
    # TotalProtein 1 g/L = 10 g/dL
    # Glucose 1 mmol/L = 18 mg/dL
    # Hemoglobin 1 g/L = 10 g/dL
    # HDL-C 1 mmol/L = 38.66976 mg/dL
    # LDL-C 1 mmol/L = 38.66976 mg/dL
    # Calcium 1 mg/dL = 0.25 mmol/L
    # UreaNitrogen 2.801 mg/dL = 1 mmol/L Urea
    for patient_id in feature_dict:
        for visit_id in feature_dict[patient_id]:
            try:
                if float(feature_dict[patient_id][visit_id]['SCr']) >= 0:
                    feature_dict[patient_id][visit_id]['SCr'] = float(feature_dict[patient_id][visit_id]['SCr']) * 88.41
            except ValueError:
                pass
            try:
                if float(feature_dict[patient_id][visit_id]['Triglycerides']) >= 0:
                    feature_dict[patient_id][visit_id]['Triglycerides'] = \
                        float(feature_dict[patient_id][visit_id]['Triglycerides']) / 88.6
            except ValueError:
                pass
            try:
                if float(feature_dict[patient_id][visit_id]['Glucose']) >= 0:
                    feature_dict[patient_id][visit_id]['Glucose'] = \
                        float(feature_dict[patient_id][visit_id]['Glucose']) / 18
            except ValueError:
                pass
            try:
                if float(feature_dict[patient_id][visit_id]['TotalBilirubin']) >= 0:
                    feature_dict[patient_id][visit_id]['TotalBilirubin'] = \
                        float(feature_dict[patient_id][visit_id]['TotalBilirubin']) * 17.1
            except ValueError:
                pass
            try:
                if float(feature_dict[patient_id][visit_id]['TotalProtein']) >= 0:
                    feature_dict[patient_id][visit_id]['TotalProtein'] = \
                        float(feature_dict[patient_id][visit_id]['TotalProtein']) * 10
            except ValueError:
                pass
            try:
                if float(feature_dict[patient_id][visit_id]['HDL-C']) >= 0:
                    feature_dict[patient_id][visit_id]['HDL-C'] = \
                        float(feature_dict[patient_id][visit_id]['HDL-C']) / 38.66976
            except ValueError:
                pass
            try:
                if float(feature_dict[patient_id][visit_id]['LDL-C']) >= 0:
                    feature_dict[patient_id][visit_id]['LDL-C'] = \
                        float(feature_dict[patient_id][visit_id]['LDL-C']) / 38.66976
            except ValueError:
                pass
            try:
                if float(feature_dict[patient_id][visit_id]['Calcium']) >= 0:
                    feature_dict[patient_id][visit_id]['Calcium'] = \
                        float(feature_dict[patient_id][visit_id]['Calcium']) / 4
            except ValueError:
                pass
            try:
                if float(feature_dict[patient_id][visit_id]['Hemoglobin']) >= 0:
                    feature_dict[patient_id][visit_id]['Hemoglobin'] = \
                        float(feature_dict[patient_id][visit_id]['Hemoglobin']) * 10
            except ValueError:
                pass
            try:
                if float(feature_dict[patient_id][visit_id]['UreaNitrogen']) >= 0:
                    feature_dict[patient_id][visit_id]['Urea'] = \
                        float(feature_dict[patient_id][visit_id]['UreaNitrogen']) / 2.801
                    feature_dict[patient_id][visit_id].pop('UreaNitrogen')
            except ValueError:
                pass
    return feature_dict


def delete_missing_critical_feature(feature_dict, critical_feature_list):
    delete_visit_set = set()
    for patient_id in feature_dict:
        for visit_id in feature_dict[patient_id]:
            for feature in critical_feature_list:
                try:
                    value = float(feature_dict[patient_id][visit_id][feature])
                    if value < 0:
                        delete_visit_set.add(patient_id + "_" + visit_id)
                except ValueError:
                    delete_visit_set.add(patient_id + "_" + visit_id)
    for item in delete_visit_set:
        patient_id, visit_id = item.split("_")
        feature_dict[patient_id].pop(visit_id)
    return feature_dict


def delete_juveniles(feature_dict, threshold=18):
    delete_visit_set = set()
    for patient_id in feature_dict:
        for visit_id in feature_dict[patient_id]:
            age = float(feature_dict[patient_id][visit_id]['age'])
            if age < threshold or age > 100:
                delete_visit_set.add(patient_id + "_" + visit_id)
    for item in delete_visit_set:
        patient_id, visit_id = item.split("_")
        feature_dict[patient_id].pop(visit_id)
    return feature_dict


def calculate_visit_count(visit_dict):
    count = 0
    for patient_id in visit_dict:
        for _ in visit_dict[patient_id]:
            count += 1
    return count


def delete_by_kidney_function(feature_dict, egfr_threshold):
    delete_visit_set = set()
    for patient_id in feature_dict:
        for visit_id in feature_dict[patient_id]:
            try:
                egfr = float(feature_dict[patient_id][visit_id]['egfr'])
                kidney_dysfunction = feature_dict[patient_id][visit_id]['CKD']
                if egfr <= egfr_threshold or kidney_dysfunction == '1' or kidney_dysfunction == 1:
                    delete_visit_set.add(patient_id + "_" + visit_id)
            except ValueError:
                print('egfr value error')
                continue
    for item in delete_visit_set:
        patient_id, visit_id = item.split("_")
        feature_dict[patient_id].pop(visit_id)
    return feature_dict


def delete_by_admission_reason(feature_dict):
    delete_visit_set = set()
    for patient_id in feature_dict:
        for visit_id in feature_dict[patient_id]:
            if len(feature_dict[patient_id][visit_id]) == 0:
                continue
            value = feature_dict[patient_id][visit_id]['HF']
            if value != '1':
                delete_visit_set.add(patient_id + "_" + visit_id)
    for item in delete_visit_set:
        patient_id, visit_id = item.split('_')
        feature_dict[patient_id].pop(visit_id)
    return feature_dict


def delete_visit_missing_too_much(feature_dict, patient_delete_missing_rate):
    # discard data whose missing rate of numerical data exceeds 0.3

    feature_type_dict = dict()
    for patient_id in feature_dict:
        for visit_id in feature_dict[patient_id]:
            for item in feature_dict[patient_id][visit_id]:
                feature_type_dict[item] = False
    for patient_id in feature_dict:
        for visit_id in feature_dict[patient_id]:
            for item in feature_dict[patient_id][visit_id]:
                value = feature_dict[patient_id][visit_id][item]
                if value != '0' and value != '1' and value != '-1':
                    feature_type_dict[item] = True
    numerical_feature_num = 0
    for item in feature_type_dict:
        if feature_type_dict[item]:
            numerical_feature_num += 1

    discard_visit_set = set()
    for patient_id in feature_dict:
        for visit_id in feature_dict[patient_id]:
            missing_count = 0
            for item in feature_dict[patient_id][visit_id]:
                if feature_type_dict[item]:
                    value = feature_dict[patient_id][visit_id][item]
                    if value == '-1':
                        missing_count += 1
            if missing_count / numerical_feature_num > patient_delete_missing_rate:
                discard_visit_set.add(patient_id + '_' + visit_id)

    for item in discard_visit_set:
        patient_id, visit_id = item.split('_')
        feature_dict[patient_id].pop(visit_id)
    return feature_dict


def read_un_preprocessed_data(source_file_path, index_item_dict):
    patient_feature_dict = dict()
    with open(source_file_path, 'r', encoding='utf-8-sig', newline="") as file:
        csv_reader = csv.reader(file)
        for line in islice(csv_reader, 2, None):
            patient_id, visit_id = line[0], line[1]
            if not patient_feature_dict.__contains__(patient_id):
                patient_feature_dict[patient_id] = dict()
            patient_feature_dict[patient_id][visit_id] = dict()
            for index in range(2, len(line)):
                item_name = index_item_dict[index]
                patient_feature_dict[patient_id][visit_id][item_name] = line[index]

    return patient_feature_dict


def discard_illegal_data_value(feature_dict):
    # 对非蛋白尿的所有数据，将数据不是数值结果的，全部置为-1
    # 对蛋白尿而言，凡出现阴性/neg，-的视为0，可数值化且数值大于5的视为1，弱阳性视为1，定性标记中存在+号的视为1，其余默认为0
    for patient_id in feature_dict:
        for visit_id in feature_dict[patient_id]:
            for item in feature_dict[patient_id][visit_id]:
                value = str(feature_dict[patient_id][visit_id][item])
                result_list = re.findall('[-+]?[\d]+(?:,\d\d\d)*[.]?\d*(?:[eE][-+]?\d+)?', value)
                if len(result_list) > 0 and result_list[0] != '-1' and result_list[0] != '-1.0':
                    feature_dict[patient_id][visit_id][item] = result_list[0]
                else:
                    feature_dict[patient_id][visit_id][item] = '-1'
    return feature_dict


def get_item_index_dict(source_file_path):
    index_item_dict = dict()
    with open(source_file_path, 'r', encoding='utf-8-sig', newline="") as file:
        csv_reader = csv.reader(file)
        for line in csv_reader:
            for index in range(2, len(line)):
                index_item_dict[index] = line[index]
            break
    return index_item_dict


if __name__ == "__main__":
    main()
