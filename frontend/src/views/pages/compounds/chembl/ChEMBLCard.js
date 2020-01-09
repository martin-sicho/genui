import { Button, CardBody, CardFooter, CardHeader, CardSubtitle } from 'reactstrap';
import React from 'react';
import './styles.css';
import TabWidget from '../../../../genui/components/TabWidget';
import ChEMBLInfo from './ChEMBLInfo';
import ChEMBLCompounds from './ChEMBLCompounds';

class ChEMBLCard extends React.Component {
  abort = new AbortController();

  constructor(props) {
    super(props);
    this.created = new Date(this.props.molset.created);
    this.updated = new Date(this.props.molset.updated);
    this.molsetURL = new URL(`chembl/${this.props.molset.id}/`, this.props.apiUrls.compoundSetsRoot);
    this.moleculesURL = new URL(`${this.props.molset.id}/molecules/`, this.props.apiUrls.compoundSetsRoot);
    this.tasksURL = new URL(`${this.props.molset.id}/tasks/all/`, this.props.apiUrls.compoundSetsRoot);

    this.state = {
      molset : null
      , isUpdating : false
      , hasUpdated : false
    }
  }

  handleDeleteSignal(deletedMolSet) {
    this.setState({isUpdating : true});
    if (this.props.hasOwnProperty("onMolsetDelete")) {
      this.props.onMolsetDelete(deletedMolSet);
    }
  }

  componentDidMount() {
    fetch(this.molsetURL, {signal : this.abort.signal})
      .then(response => this.props.handleResponseErrors(response))
      .then(this.getMolSet)
      .catch(
        (error) => console.log(error)
      )
  }

  componentWillUnmount() {
    this.abort.abort();
  }

  getMolSet = (data) => {
    this.setState({
      molset : data,
      hasUpdated : true
    })
  };

  updateMolSet = (data) => {
    this.setState({isUpdating : true});
    const error_msg = 'Failed to update ChEMBL compound set from backend.';
    fetch(
      this.molsetURL
      , {
        method: 'PATCH'
        , body: JSON.stringify(data)
        , headers: {
          'Content-Type': 'application/json'
        }
      }
    )
      .then(response => this.props.handleResponseErrors(response, error_msg))
      .then(
      data => {
        this.setState({
          molset : data,
          hasUpdated : true,
          isUpdating : false,
        })
      }
    ).catch(
      (error) => console.log(error)
    );
  };

  processTasks = (groupedTasks) => {
    if (groupedTasks.running.length > 0) {
      this.setState({
        isUpdating : true,
        hasUpdated : false
      });
    } else {
      this.setState({
        isUpdating : false,
        hasUpdated : true
      });
    }
  };

  render() {
    const molset = this.state.molset;

    if (!molset) {
      return <div>Fetching...</div>
    }

    const tabs = [
      {
        title : "Info",
        renderedComponent : () =>
          <ChEMBLInfo
            {...this.props}
            molset={molset}
            moleculesURL={this.moleculesURL}
            tasksURL={this.tasksURL}
            processTasks={this.processTasks}
            molsetIsUpdating={this.state.isUpdating}
            molsetHasUpdated={this.state.hasUpdated}
          />
      },
      {
        title: "Molecules"
        , renderedComponent : () =>
          <ChEMBLCompounds
            {...this.props}
            molset={molset}
            moleculesURL={this.moleculesURL}
            molsetIsUpdating={this.state.isUpdating}
            molsetHasUpdated={this.state.hasUpdated}
          />
      }
    ];

    return (
      <React.Fragment>
        <CardHeader>{molset.name}</CardHeader>

        <CardBody className="scrollable">
          <CardSubtitle>
            <p>
              Created: {
              this.created.toLocaleDateString()
              + ' – ' + this.created.toLocaleTimeString()
            }
              <br/>
              Last Update: {
              this.updated.toLocaleDateString()
              + ' – ' + this.updated.toLocaleTimeString()
            }
            </p>
          </CardSubtitle>
          <TabWidget tabs={tabs} />
        </CardBody>

        <CardFooter>
          <Button color="primary" disabled={this.state.isUpdating} onClick={() => this.updateMolSet({})}>{this.state.isUpdating ? 'Updating...' : 'Update Data'}</Button> <Button color="danger" disabled={this.state.isUpdating} onClick={() => {this.handleDeleteSignal(molset)}}>Delete</Button>
        </CardFooter>
      </React.Fragment>
    )
  }

}

export default ChEMBLCard;