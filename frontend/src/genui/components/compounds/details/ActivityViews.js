import React from 'react';
import { Card, CardBody, Table } from 'reactstrap';
import { groupBy, resolve } from '../../../utils';
import { TabWidget } from '../../../index';

export function ActivitiesTable(props) {
  const activities = props.activities;
  const extraData = props.extraData ? props.extraData : [];
  const extraComponents = props.extraComponents ? props.extraComponents : [];
  const appendData = props.extraDataAppend;
  const appendComponents = props.extraComponentsAppend;

  return (
    <Table size="sm" hover responsive>
      <thead>
      <tr>
        {
          !appendComponents ? extraComponents.map(component => <th key={component.header}>{component.header}</th>) : null
        }
        {
          !appendData ? extraData.map(data => <th key={data.header}>{data.header}</th>) : null
        }
        <th>Type</th>
        <th>Value</th>
        <th>Units</th>
        {
          appendComponents ? extraComponents.map(component => <th key={component.header}>{component.header}</th>) : null
        }
        {
          appendData ? extraData.map(data => <th key={data.header}>{data.header}</th>) : null
        }
      </tr>
      </thead>
      <tbody>
      {
        activities.map((activity, index) => {
          const drawExtraComponents = () => {
            return extraComponents.map(component => {
              const currentComponent = component.components[index];
              if (activity.className === currentComponent.className) {
                const Component = currentComponent.component;
                return (
                  <td key={component.header}>
                    <Component {...props} {...currentComponent.props}/>
                  </td>
                )
              } else {
                return (
                  <td key={component.header}>
                    -
                  </td>
                )
              }
            })
          };
          return (
            <tr key={activity.id}>
              {
                !appendComponents ? drawExtraComponents() : null
              }
              {
                !appendData ? extraData.map(data => <td key={data.header}>{data.data[index]}</td>) : null
              }
              <td>{activity.type.value}</td>
              <td>{activity.value.toFixed(2)}</td>
              <td>{activity.units ? activity.units.value : '-'}</td>
              {
                appendComponents ? drawExtraComponents() : null
              }
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

export function ActivitiesByTypeFlatView(props) {
  let activities = [];
  Object.keys(props.activities).forEach(key => activities = activities.concat(props.activities[key]));
  const activitySets = props.activitySets;
  if (activities.length === 0) {
    return <div>No activity data found.</div>
  }

  const activitiesGrouped = groupBy(activities, 'type.id');
  const tabs =  activitiesGrouped.map(group => {
    const extraData = [
      {
        header: "Source",
        data: group.map(activity => activitySets[activity.source].name)
      }
    ];
    const extraComponents = [];
    if (props.extraActivityFields) {
      props.extraActivityFields.forEach(definition => {
        const classNames = group.map(activity => activity.className);
        if (!classNames.includes(definition.className)) {
          return;
        }

        extraComponents.push({
          header: definition.displayName,
          components: group.map(activity => {
            // TODO: check if data items and prop names have the same length
            const data = definition.dataItems.map(dataItem => resolve(dataItem, activity));
            const sentProps = {};
            definition.propNames.forEach((propName, index) => sentProps[propName] = data[index]);
            return {
              component: definition.component,
              props: sentProps,
              className: definition.className
            }
          })
        })
      });
    }
    return {
      title: group[0].type.value,
      renderedComponent: (props) => (
        <ActivitiesTable
          {...props}
          activities={group}
          extraComponents={extraComponents}
          extraComponentsAppend={true}
          extraData={extraData}
          extraDataAppend={true}
        />
      )
    }
  });

  return (
    <TabWidget {...props} tabs={tabs}/>
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

  let sortedActivities = [];
  Object.keys(actsets).forEach(actsetKey => {
    const actset = actsets[actsetKey];
    const actsetActivities = activities[actset.id];

    if (actsetActivities &&  actsetActivities.length > 0) {
      sortedActivities = sortedActivities.concat(actsetActivities);
    }
  });
  // sortedActivities.sort((item) => item.type.name);

  const extraData = [
    {
      header: "Source",
      data: sortedActivities.map(activity => actsets[activity.source].name)
    }
  ];

  return (
    <Card>
      <CardBody>
        <ActivitiesTable
          activities={sortedActivities}
          extraData={extraData}
          extraDataAppend={true}
        />
      </CardBody>
    </Card>
  )
}