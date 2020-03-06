import React from 'react';
import { Table } from 'reactstrap';
import { groupBy } from '../../../utils';
import { TabWidget } from '../../../index';

export function ActivityTable(props) {
  const activities = props.activities;

  return (
    <Table size="sm" hover>
      <thead>
      <tr>
        <th>Value</th>
        <th>Units</th>
      </tr>
      </thead>
      <tbody>
      {
        activities.map(activity => {
          return (
            <tr key={activity.id}>
              <td>{activity.value}</td>
              <td>{activity.units ? activity.units.value : '-'}</td>
            </tr>
          )
        })
      }
      </tbody>
    </Table>
  )
}

export function ActivitiesList(props) {
  const activities = props.activities;
  const activitiesGrouped = groupBy(activities, 'type.id');
  const tabs =  activitiesGrouped.map(group => {
    return {
      title: group[0].type.value,
      renderedComponent: (props) => <ActivityTable {...props} activities={group}/>
    }
  });

  return (
    <TabWidget {...props} activitiesGrouped={activitiesGrouped} tabs={tabs}/>
  )
}

export function ActivitySetList(props) {
  const actSets = props.activitySets;
  const activities = props.activities;
  const tabs = [];
  Object.keys(actSets).forEach(key => {
    const set = actSets[key];
    if (activities[key].length > 0) {
      tabs.push({
        title: set.name,
        renderedComponent: (props) => <ActivitiesList {...props} set={set} activities={activities[key]}/>
      });
    }
  });

  return (
    <div className="activity-sets-mol-list">
      {tabs.length > 0 ? <TabWidget {...props} tabs={tabs}/> : <p>No activity information.</p>}
    </div>
  )
}