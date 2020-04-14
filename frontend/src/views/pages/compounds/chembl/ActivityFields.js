import { ComponentWithResources } from '../../../../genui';
import React from 'react';

export function AssayName(props) {
  const definition = {
    assay: new URL(`assays/${props.assayID}/`, props.molsetListUrl)
  };
  return (
    <ComponentWithResources
      {...props}
      definition={definition}
    >
      {
        (loaded, assayData) => {
          assayData = assayData.assay;
          return loaded ? <a rel="noopener noreferrer" href={`https://www.ebi.ac.uk/chembl/assay_report_card/${assayData.assayID}/`} target="_blank">{assayData.assayID}</a> : "-"
        }
      }
    </ComponentWithResources>
  )
}

export function TargetName(props) {
  const definition = {
    target: new URL(`targets/${props.targetID}/`, props.molsetListUrl)
  };
  return (
    <ComponentWithResources
      {...props}
      definition={definition}
    >
      {
        (loaded, targetData) => {
          targetData = targetData.target;
          return loaded ? <a rel="noopener noreferrer" href={`https://www.ebi.ac.uk/chembl/target_report_card/${targetData.targetID}/`} target="_blank">{targetData.targetID}</a> : "-"
        }
      }
    </ComponentWithResources>
  )
}

export function Relation(props) {
  return <span>{props.relation}</span>
}

export function ChEMBLID(props) {
  return <a rel="noopener noreferrer" href={`https://www.ebi.ac.uk/chembl/compound_report_card/${props.compoundID}/`} target="_blank">{props.compoundID}</a>
}