import csv
import os
from itertools import islice
import math


def main():
    phenogroup_path = os.path.abspath('../resource/phenogroup_assignment.csv')
    data_path = os.path.abspath('../resource/preprocessed_data.csv')
    rule_path = os.path.abspath('../resource/reproduce_mapping/value_trans_rule.csv')
    save_path = os.path.abspath('../resource/recovered_data_with_group.csv')
    phenogroup_dict = read_phenogroup_assignment(phenogroup_path)
    data_dict = read_data(data_path)
    rule_dict = read_rule_set(rule_path)
    recovered_data_dict = value_recover(data_dict, rule_dict)
    save_data(recovered_data_dict, phenogroup_dict, save_path)


def save_data(data_dict, phenogroup_dict, path):
    data_to_write = []
    head = ['patient_id', 'visit_id', 'group_id']
    for identifier in data_dict:
        for feature in data_dict[identifier]:
            head.append(feature)
        break
    data_to_write.append(head)

    for identifier in data_dict:
        line = [identifier.split('_')[0], identifier.split('_')[1], phenogroup_dict[identifier]]
        for index in range(3, len(head)):
            line.append(data_dict[identifier][head[index]])
        data_to_write.append(line)
    with open(path, 'w', encoding='utf-8-sig', newline='') as file:
        csv.writer(file).writerows(data_to_write)


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


def value_recover(data_dict, rule_dict):

    def _recover(origin_value, feature_name, recover_rule_dict):
        if not rule_dict.__contains__(feature_name):
            return origin_value
        mean = float(recover_rule_dict[feature_name]['mean'])
        stddev = float(recover_rule_dict[feature_name]['stddev'])
        value_ = origin_value * stddev + mean

        transform = recover_rule_dict[feature_name]['transform']

        try:
            if transform == 'skip':
                pass
            elif transform == 'sqrt':
                value_ = value_ ** 2
            elif transform == 'log':
                value_ = math.exp(value_)
            elif transform == 'arcsin':
                value_ = math.sin(value_ ** 2)
            else:
                raise ValueError('')
        except ZeroDivisionError:
            value_ = 0
            print('value recover error, origin value: {}, feature: {}'.format(origin_value, feature_name))

        origin_max_value = float(recover_rule_dict[feature_name]['origin_max_value'])
        origin_min_value = float(recover_rule_dict[feature_name]['origin_min_value'])

        value_ = value_ * (origin_max_value - origin_min_value) + origin_min_value
        if value_ < origin_min_value:
            value_ = origin_min_value
        if value_ > origin_max_value:
            value_ = origin_max_value

        return value_

    recovered_value_dict = dict()
    for identifier in data_dict:
        recovered_value_dict[identifier] = dict()
        for feature in data_dict[identifier]:
            value = data_dict[identifier][feature]
            recovered_value_dict[identifier][feature] = _recover(value, feature, rule_dict)
    return recovered_value_dict


def read_data(path):
    data_dict = dict()
    index_name_dict = dict()
    with open(path, 'r', encoding='utf-8-sig', newline='') as file:
        csv_reader = csv.reader(file)
        line_index = 0
        for line in csv_reader:
            if line_index == 0:
                line_index += 1
                for index in range(2, len(line)):
                    index_name_dict[index] = line[index]
                continue
            identifier = line[0] + '_' + line[1]
            data_dict[identifier] = dict()
            for index in range(2, len(line)):
                data_dict[identifier][index_name_dict[index]] = float(line[index])
    return data_dict


def read_phenogroup_assignment(path):
    phenogroup_dict = dict()
    with open(path, 'r', encoding='utf-8-sig', newline='') as file:
        csv_reader = csv.reader(file)
        for line in islice(csv_reader, 1, None):
            patient_id, visit_id, phenogroup = line[0], line[1], int(line[2])
            phenogroup_dict[patient_id+"_"+visit_id] = phenogroup
    return phenogroup_dict


if __name__ == '__main__':
    main()
