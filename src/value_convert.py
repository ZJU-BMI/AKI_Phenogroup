import csv
from itertools import islice
import os
import math


def main():
    rule_path = os.path.abspath('../resource/reproduce_mapping/value_trans_rule.csv')
    impute_path = os.path.abspath('../resource/reproduce_mapping/centroid.csv')
    data_path = os.path.abspath('../resource/filtered.csv')
    save_path = os.path.abspath('../resource/preprocessed_data.csv')

    rule_dict = read_rule_set(rule_path)
    impute_set = read_imputing_set(impute_path)
    feature_dict = read_data(data_path)
    feature_dict = value_transform(feature_dict, impute_set, rule_dict)
    data_to_write = []
    head = ['patient_id', 'visit_id']
    for patient_id in feature_dict:
        for visit_id in feature_dict[patient_id]:
            for item in feature_dict[patient_id][visit_id]:
                head.append(item)
            break
        break
    data_to_write.append(head)
    for patient_id in feature_dict:
        for visit_id in feature_dict[patient_id]:
            line = [patient_id, visit_id]
            for item in head[2:]:
                line.append(feature_dict[patient_id][visit_id][item])
            data_to_write.append(line)
    with open(save_path, 'w', encoding='utf-8-sig', newline='') as file:
        csv.writer(file).writerows(data_to_write)


def value_transform(feature_dict, impute_set, rule_dict):
    new_feature_dict = dict()
    for patient_id in feature_dict:
        new_feature_dict[patient_id] = dict()
        for visit_id in feature_dict[patient_id]:
            new_feature_dict[patient_id][visit_id] = dict()
            for item in feature_dict[patient_id][visit_id]:
                value = float(feature_dict[patient_id][visit_id][item])
                if not rule_dict.__contains__(item):
                    new_feature_dict[patient_id][visit_id][item] = value
                    continue

                if value < 0:
                    new_feature_dict[patient_id][visit_id][item] = impute_set[item]
                else:
                    max_value = float(rule_dict[item]['origin_max_value'])
                    min_value = float(rule_dict[item]['origin_min_value'])
                    transform_method = rule_dict[item]['transform']
                    if value >= max_value:
                        value = max_value - 0.000001
                    elif value <= min_value:
                        value = min_value + 0.000001
                    value = (value - min_value) / (max_value - min_value)

                    if transform_method == 'skip':
                        pass
                    elif transform_method == 'sqrt':
                        value = math.sqrt(value)
                    elif transform_method == 'log':
                        value = math.log(value)
                    elif transform_method == 'arcsin':
                        value = math.asin(value) ** 0.5
                    else:
                        raise ValueError('')

                    mean = float(rule_dict[item]['mean'])
                    stddev = float(rule_dict[item]['stddev'])
                    value = (value - mean) / stddev
                    new_feature_dict[patient_id][visit_id][item] = value
    return new_feature_dict


def read_data(source_file_path):
    patient_feature_dict = dict()
    index_item_dict = dict()
    with open(source_file_path, 'r', encoding='utf-8-sig', newline="") as file:
        csv_reader = csv.reader(file)
        line_index = 0
        for line in csv_reader:
            if line_index == 0:
                line_index += 1
                for index in range(2, len(line)):
                    index_item_dict[index] = line[index]
                continue
            patient_id, visit_id = line[0], line[1]
            if not patient_feature_dict.__contains__(patient_id):
                patient_feature_dict[patient_id] = dict()
            patient_feature_dict[patient_id][visit_id] = dict()
            for index in range(2, len(line)):
                item_name = index_item_dict[index]
                patient_feature_dict[patient_id][visit_id][item_name] = float(line[index])
    return patient_feature_dict


def read_rule_set(rule_file_path):
    recover_rule_dict = dict()
    with open(rule_file_path, 'r', encoding='utf-8-sig', newline='') as file:
        csv_reader = csv.reader(file)
        for line in islice(csv_reader, 1, None):
            recover_rule_dict[line[0]] = {
                'stddev': line[5],
                'mean': line[4],
                'origin_max_value': line[3],
                'origin_min_value': line[2],
                'transform': line[1]
            }
    return recover_rule_dict


def read_imputing_set(impute_file_path):
    impute_dict = dict()
    with open(impute_file_path, 'r', encoding='utf-8-sig', newline='') as file:
        csv_reader = csv.reader(file)
        for line in islice(csv_reader, 1, None):
            impute_dict[line[0]] = float(line[3])
    return impute_dict



if __name__ == '__main__':
    main()
