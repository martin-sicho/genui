import React from "react";
import { ComponentWithObjects, ComponentWithResources } from '../../../genui';
import ModelGrid from './ModelGrid';
import { DropdownItem, DropdownMenu, DropdownToggle, UncontrolledDropdown } from 'reactstrap';

function HeaderNav(props) {
  return (<UncontrolledDropdown nav inNavbar>
    <DropdownToggle nav caret>
      Actions
    </DropdownToggle>
    <DropdownMenu right>
      <UncontrolledDropdown>
        <DropdownToggle nav>Add New...</DropdownToggle>
        <DropdownMenu>
          {
            props.addChoices.map(choice =>
              (<DropdownItem
                key={choice.id}
                onClick={() => {props.onModelAdd(choice)}}
              >
                {choice.name}
              </DropdownItem>)
            )
          }
        </DropdownMenu>
      </UncontrolledDropdown>
    </DropdownMenu>
  </UncontrolledDropdown>)
}

class ModelsPage extends React.Component {

  constructor(props) {
    super(props);

    this.state = {
      selectedToAdd : null,
    }
  }

  handleAddNew = (model) => {
    this.setState({selectedToAdd : model})
  };

  componentDidMount() {
    this.props.onHeaderChange(<HeaderNav {...this.props} addChoices={this.props.algorithmChoices} onModelAdd={this.handleAddNew}/>);
  }

  render() {
    return (
      <div className="models-grid">
        <ComponentWithObjects
          {...this.props}
          emptyClassName="QSARModel"
          objectListURL={new URL('models/', this.props.apiUrls.qsarRoot)}
          render={
            (models, handleAddModelList, handleAddModel, handleModelDelete) => {
              return <ModelGrid
                {...this.props}
                descriptors={this.props.descriptorChoices}
                metrics={this.props.metricsChoices}
                models={models}
                chosenAlgorithm={this.state.selectedToAdd}
                handleAddModel={
                  (...args) => {
                    this.setState({selectedToAdd : null});
                    return handleAddModel(...args)
                  }
                }
                handleModelDelete={handleModelDelete}
              />
            }
          }
        />
      </div>
    );
  }
}

function Models(props) {
  const resources = {
    algorithmChoices : new URL('algorithms/', props.apiUrls.qsarRoot),
    metricsChoices: new URL('metrics/', props.apiUrls.qsarRoot),
    descriptorChoices: new URL('descriptors/', props.apiUrls.qsarRoot)
  };
  return (
    <ComponentWithResources definition={resources}>
      {
        (allLoaded, data) => (
          allLoaded ? <ComponentWithObjects
            objectListURL={new URL('all/', props.apiUrls.compoundSetsRoot)}
            {...props}
            render={
              (
                ...args
              ) => {
                const [compoundSets] = [...args];
                return (<ModelsPage
                  {...props}
                  {...data}
                  compoundSets={compoundSets}
                />)
              }
            }
          /> : <div>Loading...</div>
        )
      }
    </ComponentWithResources>
  );
}

export default Models;