import csv
import os
from itertools import islice


def main():
    phenogroup_centroid_path = os.path.abspath('../resource/reproduce_mapping/centroid.csv')
    data_path = os.path.abspath('../resource/preprocessed_data.csv')
    save_path = os.path.abspath('../resource/phenogroup_assignment.csv')
    centroid_dict = read_centroid(phenogroup_centroid_path)
    data_dict = read_data(data_path)
    assign_dict = phenogroup_assignment(centroid_dict, data_dict)
    data_to_write = [['patient_id', 'visit_id', 'phenogroup']]
    for patient_id in assign_dict:
        for visit_id in assign_dict[patient_id]:
            data_to_write.append([patient_id, visit_id, assign_dict[patient_id][visit_id]])
    with open(save_path, 'w', encoding='utf-8-sig', newline='') as file:
        csv.writer(file).writerows(data_to_write)


def phenogroup_assignment(centroid_dict, data_dict):
    phenogroup_assignment_dict = dict()
    for patient_id in data_dict:
        phenogroup_assignment_dict[patient_id] = dict()
        for visit_id in data_dict[patient_id]:
            distance_to_1 = 0
            distance_to_2 = 0
            for item in data_dict[patient_id][visit_id]:
                value = data_dict[patient_id][visit_id][item]
                distance_to_1 += (value - centroid_dict['0'][item]) ** 2
                distance_to_2 += (value - centroid_dict['1'][item]) ** 2
            if distance_to_2 < distance_to_1:
                phenogroup_assignment_dict[patient_id][visit_id] = 1
            else:
                phenogroup_assignment_dict[patient_id][visit_id] = 0
    return phenogroup_assignment_dict


def read_centroid(phenogroup_centroid_path):
    centroid_dict = {'0': dict(), '1': dict()}
    with open(phenogroup_centroid_path, 'r', encoding='utf-8-sig', newline='') as file:
        csv_reader = csv.reader(file)
        for line in islice(csv_reader, 1, None):
            feature, group1, group2 = line[0], float(line[1]), float(line[2])
            centroid_dict['0'][feature] = group1
            centroid_dict['1'][feature] = group2
    return centroid_dict


def read_data(data_path):
    patient_feature_dict = dict()
    index_item_dict = dict()
    with open(data_path, 'r', encoding='utf-8-sig', newline="") as file:
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


if __name__ == '__main__':
    main()
