import csv
import os
from joblib import load
from itertools import islice
from sklearn.metrics import roc_auc_score, recall_score
import numpy as np


def main():
    setting_dict = {
        'aki_only_group': ['group_id'],
        'aki_10':
            ['PCI', 'age', 'Acute HF', 'sex', 'antiplatelet', 'diabetes', 'NT-pro-BNP', 'vasodilator', 'Angiography',
             'CCB'],
        'aki_10_group':
            ['ACEI/ARB', 'PCI', 'egfr', 'group_id', 'Acute HF', 'sex', 'Anticoagulants', 'antiplatelet', 'TnT',
             'NT-pro-BNP', 'Angiography'],
        'aki_39':
            ['ACEI/ARB', 'BMI', 'PCI', 'beta-blocker', 'egfr', 'GGT', 'ALT', 'LDL-C', 'CHD', 'diuretic', 'AST',
             'LVEF', 'Urea', 'age', 'AF', 'myocardiopathy', 'Acute HF', 'sex', 'TotalBilirubin', 'TotalProtein',
             'Anticoagulants', 'antiplatelet', 'PositiveInotropicDrugs', 'VHD', 'Triglycerides', 'diabetes', 'TnT',
             'NT-pro-BNP', 'stroke', 'DBP', 'SBP', 'vasodilator', 'Hemoglobin', 'Angiography', 'Calcium', 'CCB',
             'Sodium', 'Potassium', 'HDL-C'],
        'aki_39_group':
            ['ACEI/ARB', 'BMI', 'PCI', 'beta-blocker', 'egfr', 'group_id', 'GGT', 'ALT', 'LDL-C', 'CHD', 'diuretic',
             'AST', 'LVEF', 'Urea', 'age', 'AF', 'myocardiopathy', 'Acute HF', 'sex', 'TotalBilirubin', 'TotalProtein',
             'Anticoagulants', 'antiplatelet', 'PositiveInotropicDrugs', 'VHD', 'Triglycerides', 'diabetes', 'TnT',
             'NT-pro-BNP', 'stroke', 'DBP', 'SBP', 'vasodilator', 'Hemoglobin', 'Angiography', 'Calcium', 'CCB',
             'Sodium', 'Potassium', 'HDL-C'],
        'death_only_group': ['group_id'],
        'death_10':
            ['PCI', 'egfr', 'AST', 'age', 'Acute HF', 'sex', 'TotalBilirubin', 'diabetes', 'vasodilator',
             'Angiography'],
        'death_10_group':
            ['PCI', 'egfr', 'group_id', 'ALT', 'AST', 'Urea', 'age', 'Acute HF', 'sex', 'diabetes', 'Angiography'],
        'death_39':
            ['ACEI/ARB', 'BMI', 'PCI', 'beta-blocker', 'egfr', 'GGT', 'ALT', 'LDL-C', 'CHD', 'diuretic', 'AST',
             'LVEF', 'Urea', 'age', 'AF', 'myocardiopathy', 'Acute HF', 'sex', 'TotalBilirubin', 'TotalProtein',
             'Anticoagulants', 'antiplatelet', 'PositiveInotropicDrugs', 'VHD', 'Triglycerides', 'diabetes', 'TnT',
             'NT-pro-BNP', 'stroke', 'DBP', 'SBP', 'vasodilator', 'Hemoglobin', 'Angiography', 'Calcium', 'CCB',
             'Sodium', 'Potassium', 'HDL-C'],
        'death_39_group':
            ['ACEI/ARB', 'BMI', 'PCI', 'beta-blocker', 'egfr', 'group_id', 'GGT', 'ALT', 'LDL-C', 'CHD', 'diuretic',
             'AST', 'LVEF', 'Urea', 'age', 'AF', 'myocardiopathy', 'Acute HF', 'sex', 'TotalBilirubin', 'TotalProtein',
             'Anticoagulants', 'antiplatelet', 'PositiveInotropicDrugs', 'VHD', 'Triglycerides', 'diabetes', 'TnT',
             'NT-pro-BNP', 'stroke', 'DBP', 'SBP', 'vasodilator', 'Hemoglobin', 'Angiography', 'Calcium', 'CCB',
             'Sodium', 'Potassium', 'HDL-C'],
    }

    label_path = os.path.abspath('../resource/mimic_unpreprocessed.csv')
    data_path = os.path.abspath('../resource/recovered_data_with_group.csv')
    model_folder = os.path.abspath('../resource/reproduce_model')
    aki_dict, death_dict = read_label(label_path, aki_index=9, death_index=13)
    data_dict = read_origin_data(data_path)

    for experiment_name in setting_dict:
        if experiment_name.__contains__('death'):
            label_dict = death_dict
        elif experiment_name.__contains__('aki'):
            label_dict = aki_dict
        else:
            raise ValueError('')
        prediction_experiment(label_dict, setting_dict[experiment_name], experiment_name, model_folder, data_dict)


def prediction_experiment(label_dict, feature_list, model_name, model_folder, data_dict):
    predict_model = load(os.path.join(model_folder, model_name+'.joblib'))
    data_list, label_list = reconstruct_data(data_dict, feature_list, label_dict)
    predict_label = predict_model.predict(data_list)
    predict_prob = predict_model.predict_proba(data_list)[:, 1]

    auc = roc_auc_score(label_list, predict_prob)
    recall = recall_score(label_list, predict_label)
    specificity = np.sum((-1 * label_list + 1) * (-1 * predict_label + 1)) / (
            len(label_list) - np.sum(label_list))
    print('{}, recall (sensitivity): {:.5f}, specificity: {:.5f}, auc: {:.5f}'
          .format(model_name, recall, specificity, auc))


def reconstruct_data(data_dict, feature_list, label_dict):
    data_list = []
    label_list = []
    for identifier in data_dict:
        line = []
        for feature in feature_list:
            line.append(data_dict[identifier][feature])
        data_list.append(line)
        label_list.append(label_dict[identifier])
    return np.array(data_list), np.array(label_list)


def read_origin_data(path):
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


def read_label(path, aki_index, death_index):
    aki_dict = dict()
    death_dict = dict()
    with open(path, 'r', encoding='utf-8-sig', newline='') as file:
        csv_reader = csv.reader(file)
        for line in islice(csv_reader, 1, None):
            identifier = line[0]+"_"+line[1]
            aki_dict[identifier] = int(line[aki_index])
            death_dict[identifier] = int(line[death_index])
    return aki_dict, death_dict


def read_model():
    pass


if __name__ == '__main__':
    main()
