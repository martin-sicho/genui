import React from 'react';
import { GenericMolSetCard, MolsetActivitiesSummary, MolsInMolSetList } from '../../../../genui';
import GeneratedInfo from './tabs/GeneratedInfo';

function GeneratedCard(props) {
  const tabs = [
    {
      title : "Info",
      renderedComponent : GeneratedInfo,
    },
    {
      title : "Structures",
      renderedComponent: (props) => <MolsInMolSetList {...props} showInfo={true}/>,
    },
    {
      title: "Activities",
      renderedComponent: props => <MolsetActivitiesSummary {...props} selectable={false}/>
    }
  ];

  return (
    <GenericMolSetCard {...props} tabs={tabs}/>
  )
}

export default GeneratedCard;