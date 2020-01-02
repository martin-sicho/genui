import { Button, CardBody, CardFooter, CardHeader, CardSubtitle } from 'reactstrap';
import React from 'react';
import './styles.css';
import TabWidget from './TabWidget';
import ChEMBLInfo from './ChEMBLInfo';
import ChEMBLCompounds from './ChEMBLCompounds';

class ChEMBLCard extends React.Component {

  constructor(props) {
    super(props);
    this.created = new Date(this.props.molset.created);
    this.updated = new Date(this.props.molset.updated);
    this.molsetURL = new URL(`chembl/${this.props.molset.id}/`, this.props.apiUrls.compoundSetsRoot);
    this.moleculesURL = new URL(`${this.props.molset.id}/molecules/`, this.props.apiUrls.compoundSetsRoot);
    this.tasksURL = new URL(`${this.props.molset.id}/tasks/all/`, this.props.apiUrls.compoundSetsRoot);

    this.state = {
      molset : null
    }
  }

  componentDidMount() {
    fetch(this.molsetURL)
      .then(response => response.json())
      .then(this.getMolSet)
  }

  getMolSet = (data) => {
    this.setState({molset : data})
  };

  render() {
    const molset = this.state.molset;

    if (!molset) {
      return <div>Fetching...</div>
    }

    const tabs = [
      {
        title : "Info",
        renderedComponent : () => <ChEMBLInfo {...this.props} molset={molset} moleculesURL={this.moleculesURL} tasksURL={this.tasksURL} />
      },
      {
        title: "Molecules"
        , renderedComponent : () => <ChEMBLCompounds {...this.props} molset={molset} moleculesURL={this.moleculesURL} />
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

        {/*TODO: implement*/}
        <CardFooter>
          <Button color="primary">Update Data</Button> <Button color="danger">Delete</Button>
        </CardFooter>
      </React.Fragment>
    )
  }

}

export default ChEMBLCard;