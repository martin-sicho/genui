import React from 'react';
import { ChEMBLCreateForm } from './CreateForm';
import { CardHeader } from 'reactstrap';

class ChEMBLCardNew extends React.Component {

  createMolSetFromFormData = (data) => {
    data.project = this.props.currentProject.id;
    if (!data.maxPerTarget) delete data.maxPerTarget;
    const url = new URL('chembl/', this.props.apiUrls.compoundSetsRoot);
    fetch(
      url
      , {
        method: 'POST'
        , body: JSON.stringify(data)
        , headers: {
          'Content-Type': 'application/json'
        }
      }
    ).then(response => response.json()).then(
      data => {
        this.props.handleCreateNew(this.props.currentMolsetClass, data)
      }
    );
  };

  render() {
    return (
      <React.Fragment>
        <CardHeader>Download New Compound Set from ChEMBL</CardHeader>
        <ChEMBLCreateForm handleCreate={this.createMolSetFromFormData}/>
      </React.Fragment>
    )
  }

}

export default ChEMBLCardNew;