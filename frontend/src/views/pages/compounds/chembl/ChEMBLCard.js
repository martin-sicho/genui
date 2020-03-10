import React from 'react';
import ChEMBLInfo from './tabs/ChEMBLInfo';
import { GenericMolSetCard, MolsInMolSetList } from '../../../../genui';

function ChEMBLCard(props) {
  const tabs = [
    {
      title : "Info",
      renderedComponent : ChEMBLInfo,
    },
    {
      title: "Molecules",
      renderedComponent : (props) => <MolsInMolSetList {...props} showInfo={false}/>
    }
  ];

  return (
    <GenericMolSetCard {...props} tabs={tabs}/>
  )
}

export default ChEMBLCard;