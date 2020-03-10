import React from 'react';
import { GenericMolSetCard, MolsInMolSetList } from '../../../../genui';
import GeneratedInfo from './tabs/GeneratedInfo';

function GeneratedCard(props) {
  const tabs = [
    {
      title : "Info",
      renderedComponent : GeneratedInfo,
    },
    {
      title : "Molecules",
      renderedComponent: (props) => <MolsInMolSetList {...props} showInfo={false}/>,
    }
  ];

  return (
    <GenericMolSetCard {...props} tabs={tabs}/>
  )
}

export default GeneratedCard;