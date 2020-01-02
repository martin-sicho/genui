import { Button, CardBody, CardFooter, CardHeader, CardSubtitle } from 'reactstrap';
import React from 'react';
import './styles.css';
import TabWidget from './TabWidget';
import ChEMBLInfo from './ChEMBLInfo';
import ChEMBLCompounds from './ChEMBLCompounds';

class ChEMBLCard extends React.Component {

  constructor(props) {
    super(props);
    this.molset = props.molset;
    this.created = new Date(this.molset.created);
    this.updated = new Date(this.molset.updated);
  }

  render() {
    const tabs = [
      {
        title : "Info",
        renderedComponent : () => <ChEMBLInfo molset={this.molset}/>
      },
      {
        title: "Molecules"
        , renderedComponent : () => <ChEMBLCompounds molset={this.molset}/>
      }
    ];

    return (
      <React.Fragment>
        <CardHeader>{this.molset.name}</CardHeader>
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
          <Button color="success">Add</Button> <Button>Cancel</Button>
        </CardFooter>
      </React.Fragment>
    )
  }

}

export default ChEMBLCard;