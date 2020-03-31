import React from 'react';
import { Card, CardBody, CardHeader, Table } from 'reactstrap';
import { groupBy } from '../../../utils';
import { TabWidget } from '../../../index';

export function ActivitiesTable(props) {
  const activities = props.activities;
  const extraData = props.extraData ? props.extraData : [];
  const appendData = props.extraDataAppend;

  return (
    <Table size="sm" hover>
      <thead>
      <tr>
        {
          !appendData ? extraData.map(data => <th key={data.header}>{data.header}</th>) : null
        }
        <th>Value</th>
        <th>Units</th>
        {
          appendData ? extraData.map(data => <th key={data.header}>{data.header}</th>) : null
        }
      </tr>
      </thead>
      <tbody>
      {
        activities.map((activity, index) => {
          return (
            <tr key={activity.id}>
              {
                !appendData ? extraData.map(data => <td key={data.header}>{data.data[index]}</td>) : null
              }
              <td>{activity.value.toFixed(2)}</td>
              <td>{activity.units ? activity.units.value : '-'}</td>
              {
                appendData ? extraData.map(data => <td key={data.header}>{data.data[index]}</td>) : null
              }
            </tr>
          )
        })
      }
      </tbody>
    </Table>
  )
}

export function ActivitiesByTypeTabView(props) {
  const activities = props.activities;
  const activitiesGrouped = groupBy(activities, 'type.id');
  const tabs =  activitiesGrouped.map(group => {
    return {
      title: group[0].type.value,
      renderedComponent: (props) => <ActivitiesTable {...props} activities={group}/>
    }
  });

  return (
    <TabWidget {...props} activitiesGrouped={activitiesGrouped} tabs={tabs}/>
  )
}

export function ActivitySetTabView(props) {
  const actSets = props.activitySets;
  const activities = props.activities;
  const tabs = [];
  Object.keys(actSets).forEach(key => {
    const set = actSets[key];
    if (activities[key].length > 0) {
      tabs.push({
        title: set.name,
        renderedComponent: (props) => <ActivitiesByTypeTabView {...props} set={set} activities={activities[key]}/>
      });
    }
  });

  return (
    <div className="activity-sets-mol-list">
      {tabs.length > 0 ? <TabWidget {...props} tabs={tabs}/> : <p>No activity information.</p>}
    </div>
  )
}

export function ActivitySetFlatView(props) {
  const activities = props.activities;
  const actsets = props.activitySets;

  return (
    <React.Fragment>
      {
        Object.keys(actsets).map(actsetKey => {
          const actset = actsets[actsetKey];
          const actsetActivities = activities[actset.id];
          if (actsetActivities &&  actsetActivities.length > 0) {
            // const byType = groupBy(actsetActivities, 'type.id');
            // console.log(byType);
            actsetActivities.sort((item) => item.type.name);

            // FIXME: filtering should not be necessary, fix ComponentWithPagedResources or MoleculeActivityProvider so that it does not leak previous information
            const filteredActivities = [];
            actsetActivities.forEach(activity => {
              if (activity.molecule === props.mol.id && !filteredActivities.find(item => item.id === activity.id)) {
                filteredActivities.push(activity);
              }
            });
            const extraData = [
              {
                header: "Type",
                data: filteredActivities.map(activity => activity.type.value)
              }
            ];
            return (
              <Card key={actsetKey}>
                <CardHeader>{actset.name}</CardHeader>
                <CardBody>
                  <ActivitiesTable
                    activities={filteredActivities}
                    extraData={extraData}
                    extraDataAppend={false}
                  />
                </CardBody>
              </Card>
            )
          } else {
            return null;
          }
        })
      }
    </React.Fragment>
  )
}