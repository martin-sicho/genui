import React from "react";
import {
  ComponentWithObjects,
  ComponentWithResources,
  ModelsPage,
} from '../../../genui';
import QSARModelCard from './ModelCard';
import QSARModelCreateCard from './QSARModelCreateCard';
import { DropdownItem, DropdownMenu, DropdownToggle, UncontrolledDropdown } from 'reactstrap';
import QSARModelCreateFromFileCard from './QSARModelCreateFromFileCard';

function HeaderNav(props) {
  return (<UncontrolledDropdown nav inNavbar>
    <DropdownToggle nav caret>
      Create QSAR Model
    </DropdownToggle>
    <DropdownMenu right>
      <UncontrolledDropdown>
        <DropdownToggle nav>Train...</DropdownToggle>
        <DropdownMenu>
          {
            props.addChoices.map(choice =>
              (<DropdownItem
                key={choice.id}
                onClick={() => {props.onModelAdd(choice, QSARModelCreateCard, {
                  h : {"md" : 10, "sm" : 7},
                  w : {"md" : 1, "sm" : 1},
                  minH : {"md" : 3, "sm" : 3},
                })}}
              >
                {choice.name}
              </DropdownItem>)
            )
          }
        </DropdownMenu>
      </UncontrolledDropdown>
      <UncontrolledDropdown>
        <DropdownToggle nav>From File...</DropdownToggle>
        <DropdownMenu>
          {
            props.addChoices.map(choice =>
              (<DropdownItem
                key={choice.id}
                onClick={() => {props.onModelAdd(choice, QSARModelCreateFromFileCard, {
                  h : {"md" : 8, "sm" : 8},
                  w : {"md" : 1, "sm" : 1},
                  minH : {"md" : 3, "sm" : 3},
                })}}
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

function Models(props) {
  const resources = {
    algorithmChoices : new URL('algorithms/', props.apiUrls.qsarRoot),
    descriptors: new URL('descriptors/', props.apiUrls.qsarRoot),
    metrics: new URL('metrics/', props.apiUrls.qsarRoot)
  };
  const listUrl = new URL(`models/`, props.apiUrls.qsarRoot);
  const defaultClassName = "QSARModel";
  return (
    <ComponentWithResources definition={resources}>
      {
        (allLoaded, resources) => (
          allLoaded ? <ComponentWithObjects
            {...props}
            objectListURL={new URL('all/', props.apiUrls.compoundSetsRoot)}
            emptyClassName={defaultClassName}
            render={
              (
                ...args
              ) => {
                const [compoundSets] = [...args];
                const compoundSetsAvailable = !(Object.keys(compoundSets).length === 1 && compoundSets[defaultClassName].length === 0 && compoundSets.constructor === Object);
                return (compoundSetsAvailable ? <ModelsPage
                  {...props}
                  {...resources}
                  modelClass={defaultClassName}
                  listURL={listUrl}
                  modelComponent={QSARModelCard}
                  compoundSets={compoundSets}
                  headerComponent={HeaderNav}
                  cardSetup={{
                    h : {"md" : 10, "sm" : 7},
                    w : {"md" : 1, "sm" : 1},
                    minH : {"md" : 3, "sm" : 3},
                  }}
                /> : <div><p>There are currently no compound sets. You need to create one before building a QSAR model.</p></div>)
              }
            }
          /> : <div>Loading...</div>
        )
      }
    </ComponentWithResources>
  );
}

export default Models;