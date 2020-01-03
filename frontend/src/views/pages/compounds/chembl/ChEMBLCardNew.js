import React from 'react';
import { ChEMBLCreateForm } from './CreateForm';
import { CardHeader } from 'reactstrap';

class ChEMBLCardNew extends React.Component {

  createMolSetFromFormData = (data) => {
    const molset_new = data;
    this.props.handleCreateNew(molset_new)
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