from lifelines import KaplanMeierFitter
import os
import csv
from itertools import islice
import pandas as pd
import matplotlib.pyplot as plt
from lifelines.statistics import multivariate_logrank_test


def read_group_id(group_id_path):
    group_id_dict = dict()
    with open(group_id_path, 'r', encoding='utf-8-sig', newline='') as file:
        csv_reader = csv.reader(file)
        for line in islice(csv_reader, 1, None):
            patient_id, visit_id, group_id = line
            if not group_id_dict.__contains__(patient_id):
                group_id_dict[patient_id] = dict()
            group_id_dict[patient_id][visit_id] = group_id
    return group_id_dict


def read_event(file_path, aki_index, death_index, aki_time_index, death_time_index, threshold=30):
    event_dict = dict()
    with open(file_path, 'r', encoding='utf-8-sig', newline='') as file:
        csv_reader = csv.reader(file)
        for line in islice(csv_reader, 2, None):
            patient_id, visit_id, aki, aki_time, death, death_time = \
                line[0], line[1], line[aki_index], float(line[aki_time_index]), line[death_index], \
                float(line[death_time_index])
            if float(line[death_index]) > 0.5:
                death = 1
            else:
                death = 0
            if float(line[aki_index]) > 0.5:
                aki = 1
            else:
                aki = 0

            if not event_dict.__contains__(patient_id):
                event_dict[patient_id] = dict()

            if aki == 0:
                aki_time = threshold
            if aki == 1 and aki_time > threshold:
                aki_time = threshold
            if death == 0:
                death_time = threshold
            if death == 1 and death_time > threshold:
                death_time = threshold
            event_dict[patient_id][visit_id] = {'aki_event': aki, 'aki_time': aki_time, 'death_event': death,
                                                'death_time': death_time}

    return event_dict


def data_fusion(event_dict, group_id_dict):
    aki_dict = {}
    aki_time_dict = {}
    death_dict = {}
    death_time_dict = {}
    for patient_id in group_id_dict:
        for visit_id in group_id_dict[patient_id]:
            group_id = group_id_dict[patient_id][visit_id]
            if not aki_dict.__contains__(group_id):
                aki_dict[group_id] = []
                aki_time_dict[group_id] = []
                death_dict[group_id] = []
                death_time_dict[group_id] = []
            aki_event = event_dict[patient_id][visit_id]['aki_event'] == 1
            aki_time = event_dict[patient_id][visit_id]['aki_time']
            death_event = event_dict[patient_id][visit_id]['death_event']
            death_time = event_dict[patient_id][visit_id]['death_time']
            aki_dict[group_id].append(aki_event)
            aki_time_dict[group_id].append(aki_time)
            death_dict[group_id].append(death_event)
            death_time_dict[group_id].append(death_time)
    return aki_dict, aki_time_dict, death_dict, death_time_dict


def log_rank_test(event_dict, time_dict, label):
    group = []
    duration = []
    event = []
    for group_id in event_dict:
        for index in range(len(event_dict[group_id])):
            group.append(int(group_id))
            duration.append(time_dict[group_id][index])
            event.append(event_dict[group_id][index])
    df = pd.DataFrame({'durations': duration, 'groups': group, 'events': event,})
    print(label)
    results = multivariate_logrank_test(df['durations'], df['groups'], df['events'])
    results.print_summary()


def main():
    mimic_group_id_path = os.path.abspath('../resource/phenogroup_assignment.csv')
    mimic_event_file_path = os.path.abspath('../resource/mimic_unpreprocessed.csv')

    mimic_event_dict = read_event(mimic_event_file_path, aki_index=9, death_index=13, aki_time_index=10,
                                  death_time_index=14, threshold=30)
    mimic_group_id_dict = read_group_id(mimic_group_id_path)

    m_aki_dict, m_aki_time_dict, m_death_dict, m_death_time_dict = data_fusion(mimic_event_dict, mimic_group_id_dict)

    # log rank test
    log_rank_test(m_aki_dict, m_aki_time_dict, 'mimic aki')
    log_rank_test(m_death_dict, m_death_time_dict, 'mimic death')

    kmf = KaplanMeierFitter()

    # Fit the data into the model
    fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(8, 4))

    group_list = ['0', '1']
    color_list = ['#070707', '#ff3b3b']
    for group_id in group_list:
        kmf.fit(m_aki_time_dict[group_id], event_observed=m_aki_dict[group_id], label='Phenogroup {}'
                .format(int(group_id)+1)).plot(ax=axs[0], ci_show=False, color=color_list[int(group_id)],
                                               legend=False)
    for group_id in group_list:
        kmf.fit(m_death_time_dict[group_id], event_observed=m_death_dict[group_id], label='Phenogroup {}'
                .format(int(group_id)+1)).plot(ax=axs[1], ci_show=False, color=color_list[int(group_id)],
                                               legend=False)

    axs[0].text(1, 0.5416, "log rank p < 0.005", ha='left')
    axs[1].text(1, 0.89, "log rank p = 0.01", ha='left')
    axs[0].title.set_text('AKI in MIMIC Dataset')
    axs[1].title.set_text('In-Hospital Mortality in MIMIC Dataset')
    axs[0].set_xlabel('Time (days)')
    axs[1].set_xlabel('Time (days)')
    axs[0].set_ylabel('Survival free of AKI (%)')
    axs[1].set_ylabel('Survival free of Death (%)')
    axs[0].set_ylim(0.5, 1.01)
    axs[1].set_ylim(0.88, 1.005)
    axs[0].set_xlim(0, 31)
    axs[1].set_xlim(0, 31)
    handles, labels = axs[0].get_legend_handles_labels()
    legend = fig.legend(handles, labels, bbox_to_anchor=(0.74, 1.08), ncol=2)
    plt.tight_layout()

    plt.show()
    fig.savefig(os.path.abspath('../resource/kaplan_meier_curve.png'),
                bbox_extra_artists=(legend,),
                bbox_inches='tight'
                )


if __name__ == '__main__':
    main()
