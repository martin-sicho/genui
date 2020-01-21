import React from "react"
import { CardHeader } from 'reactstrap';
import ModelCreateForm from './CreateForm';

class ModelCardNew extends React.Component {

  newModelFromFormData = (data) => {
    const url = new URL('models/', this.props.apiUrls.qsarRoot);
    fetch(
      url
      , {
        method: 'POST'
        , body: JSON.stringify(data)
        , headers: {
          'Content-Type': 'application/json'
        }
      }
    ).then((data) => this.props.handleResponseErrors(data, "Bad Request. Data wrong or incomplete?")).then(
      data => {
        this.props.handleAddModel('QSARModel', data)
      }
    ).catch(
      error => console.log(error)
    );
  };

  render() {
    let molsets = [];
    Object.keys(this.props.compoundSets).forEach(
      (key) => molsets = molsets.concat(this.props.compoundSets[key])
    );

    if (molsets.length === 0) {
      // TODO: do something smarter like display an alert or redirect to the compounds UI
      return <p>There are no existing compound sets. Create one first.</p>
    }

    return (
      <React.Fragment>
        <CardHeader>Create New {this.props.chosenAlgorithm.name} Model</CardHeader>
        <ModelCreateForm
          handleCreate={this.newModelFromFormData}
          chosenAlgorithm={this.props.chosenAlgorithm}
          molsets={molsets}
          descriptors={this.props.descriptors}
          metrics={this.props.metrics}
          project={this.props.currentProject}
        />
      </React.Fragment>
    )
  }
}

export default ModelCardNew;