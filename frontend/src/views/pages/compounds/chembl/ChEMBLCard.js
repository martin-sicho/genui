import React from 'react';
import ChEMBLInfo from './tabs/ChEMBLInfo';
import { GenericMolSetCard, MolsetActivitiesSummary, MolsInMolSetList } from '../../../../genui';
import { AssayName, ChEMBLID, Relation, TargetName } from './ActivityFields';

function ChEMBLCard(props) {
  const tabs = [
    {
      title : "Info",
      renderedComponent : ChEMBLInfo,
    },
    {
      title: "Compounds",
      renderedComponent : (props) => (
        <MolsInMolSetList
          {...props}
          showInfo={true}
          extraActivityFields={
            [
              {
                dataItems: ["extraArgs.relation"],
                propNames: ["relation"],
                displayName: "Relation",
                component: Relation,
                className: "ChEMBLActivity"
              },
              {
                dataItems: ["extraArgs.assay"],
                propNames: ["assayID"],
                displayName: "Assay",
                component: AssayName,
                className: "ChEMBLActivity"
              },
              {
                dataItems: ["extraArgs.target"],
                propNames: ["targetID"],
                displayName: "Target",
                component: TargetName,
                className: "ChEMBLActivity"
              },
            ]
          }
          extraInfoFields={
            [
              {
                dataItems: ["extraArgs.chemblID"],
                propNames: ["compoundID"],
                displayName: "ChEMBL ID",
                component: ChEMBLID,
                className: "ChEMBLMolecule"
              }
            ]
          }
        />)
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

export default ChEMBLCard;