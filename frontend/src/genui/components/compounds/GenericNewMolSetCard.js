import React from "react";
import { CardHeader } from 'reactstrap';
import NewMolSetFormRenderer from './NewMolSetFormRenderer';

class GenericNewMolSetCard extends React.Component {

  createMolSetFromFormData = (data) => {
    data.project = this.props.currentProject.id;
    if (!data.maxPerTarget) delete data.maxPerTarget; // TODO: remove this and make sure nothing breaks
    fetch(
      this.props.molsetListUrl
      , {
        method: 'POST'
        , body: JSON.stringify(data)
        , headers: {
          'Content-Type': 'application/json'
        },
        credentials: "include",
      }
    ).then(response => response.json()).then(
      data => {
        this.props.handleCreateNew(this.props.currentMolsetClass, data)
      }
    ).catch(
      (e) => console.log(e)
    );
  };

  render() {
    return (
      <React.Fragment>
        <CardHeader>{this.props.cardHeader}</CardHeader>
        <NewMolSetFormRenderer {...this.props} handleCreate={this.createMolSetFromFormData}/>
      </React.Fragment>
    )
  }
}

export default GenericNewMolSetCard;